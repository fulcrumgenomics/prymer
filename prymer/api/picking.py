"""
# Methods and class for building and scoring primer pairs from a list of left and right primers.

Typically, primer pairs are built from a list of left and right primers for a
given target with the [`build_primer_pairs()`][prymer.api.picking.build_primer_pairs]
method. The returned primer pairs are automatically scored (using the
[`score()`][prymer.api.picking.score] method).

## Module contents

Contains the following public classes and methods:

- [`score()`][prymer.api.picking.score] -- Scores the amplicon amplified by a primer pair
    in a manner similar to Primer3.
- [`build_primer_pairs()`][prymer.api.picking.build_primer_pairs] -- Builds primer pairs
    from individual left and primers.

"""

from collections.abc import Sequence
from pathlib import Path
from typing import Iterator
from typing import Optional

from pysam import FastaFile

from prymer.api.melting import calculate_long_seq_tm
from prymer.api.minoptmax import MinOptMax
from prymer.api.oligo import Oligo
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.ntthal import NtThermoAlign
from prymer.primer3 import PrimerAndAmpliconWeights


def score(
    left_primer: Oligo,
    right_primer: Oligo,
    amplicon: Span,
    amplicon_tm: float,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
    weights: PrimerAndAmpliconWeights,
) -> float:
    """Score the amplicon in a manner similar to Primer3

    Sums the following:

    1. Primer penalties: the provided left and right primer penalties
    2. Amplicon size: The difference between the current amplicon size and optimal amplicon size
       scaled by the product size weight.  Is zero if the optimal amplicon size is zero.
    3. Amplicon Tm: The difference in melting temperature between the calculated and optimal,
       weighted by the product melting temperature.

    Args:
        left_primer: the left primer
        right_primer: the right primer
        amplicon: the amplicon mapping
        amplicon_tm: the melting temperature of the amplicon
        amplicon_sizes: minimum, optimal, and maximum amplicon sizes (lengths)
        amplicon_tms: minimum, optimal, and maximum amplicon Tms
        weights: the set of penalty weights

    Returns:
        the penalty for the whole amplicon.
    """
    # The penalty for the amplicon size:
    # 1. No penalty if the optimal amplicon size is zero
    # 2. The difference between the current amplicon size and optimal amplicon size scaled by the
    #    product size weight.  The product size weight is different depending on if the amplicon
    #    size is greater or less than the optimal amplicons.
    size_penalty: float
    if amplicon_sizes.opt == 0:
        size_penalty = 0.0
    elif amplicon.length > amplicon_sizes.opt:
        size_penalty = (amplicon.length - amplicon_sizes.opt) * weights.product_size_gt
    else:
        size_penalty = (amplicon_sizes.opt - amplicon.length) * weights.product_size_lt

    # The penalty for the amplicon melting temperature.
    # The difference in melting temperature between the calculated and optimal is weighted by the
    # product melting temperature.
    tm = amplicon_tm
    tm_penalty: float
    if tm > amplicon_tms.opt:
        tm_penalty = (tm - amplicon_tms.opt) * weights.product_tm_gt
    else:
        tm_penalty = (amplicon_tms.opt - tm) * weights.product_tm_lt

    # Put it all together
    return left_primer.penalty + right_primer.penalty + size_penalty + tm_penalty


def build_primer_pairs(
    left_primers: Sequence[Oligo],
    right_primers: Sequence[Oligo],
    target: Span,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
    max_heterodimer_tm: Optional[float],
    weights: PrimerAndAmpliconWeights,
    fasta_path: Path,
) -> Iterator[PrimerPair]:
    """Builds primer pairs from individual left and primers.

    Only primer pairs that meet the requirements for amplicon sizes and tms will be returned.
    In addition, if `max_heterodimer_tm` is provided then primer pairs with a heterodimer tm
    exceeding the maximum will also be discarded.

    Args:
        left_primers: the left primers
        right_primers: the right primers
        target: the genome mapping for the target
        amplicon_sizes: minimum, optimal, and maximum amplicon sizes (lengths)
        amplicon_tms: minimum, optimal, and maximum amplicon Tms
        max_heterodimer_tm: if supplied, heterodimer Tms will be calculated for primer pairs,
            and those exceeding the maximum Tm will be discarded
        weights: the set of penalty weights
        fasta_path: the path to the FASTA file from which the amplicon sequence will be retrieved.

    Returns:
        an iterator over all the valid primer pairs, unsorted
    """
    # Short circuit if we have no left primers or no right primers
    if not any(left_primers) or not any(right_primers):
        return

    if any(p.span.refname != target.refname for p in left_primers):
        raise ValueError("Left primers exist on different reference than target.")

    if any(p.span.refname != target.refname for p in right_primers):
        raise ValueError("Right primers exist on different reference than target.")

    # Grab the sequence we'll use to fill in the amplicon sequence
    with FastaFile(f"{fasta_path}") as fasta:
        region_start = min(p.span.start for p in left_primers)
        region_end = max(p.span.end for p in right_primers)
        bases = fasta.fetch(target.refname, region_start - 1, region_end)

    with NtThermoAlign() as ntthal:
        # generate all the primer pairs that don't violate hard size and Tm constraints
        for lp in left_primers:
            for rp in right_primers:
                if rp.span.end - lp.span.start + 1 > amplicon_sizes.max:
                    continue

                amp_mapping = Span(refname=target.refname, start=lp.span.start, end=rp.span.end)
                amp_bases = bases[
                    amp_mapping.start - region_start : amp_mapping.end - region_start + 1
                ]
                amp_tm = calculate_long_seq_tm(amp_bases)

                if amp_tm < amplicon_tms.min or amp_tm > amplicon_tms.max:
                    continue

                if max_heterodimer_tm is not None:
                    if ntthal.duplex_tm(lp.bases, rp.bases) > max_heterodimer_tm:
                        continue

                penalty = score(
                    left_primer=lp,
                    right_primer=rp,
                    amplicon=amp_mapping,
                    amplicon_tm=amp_tm,
                    amplicon_sizes=amplicon_sizes,
                    amplicon_tms=amplicon_tms,
                    weights=weights,
                )

                pp = PrimerPair(
                    left_primer=lp,
                    right_primer=rp,
                    amplicon_sequence=amp_bases,
                    amplicon_tm=amp_tm,
                    penalty=penalty,
                )

                yield pp
