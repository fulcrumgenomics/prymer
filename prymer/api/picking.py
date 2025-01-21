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
from typing import Tuple

from pysam import FastaFile

from prymer import Thermo
from prymer.model import MinOptMax
from prymer.model import Oligo
from prymer.model import PrimerPair
from prymer.model import Span
from prymer.primer3 import PrimerParameters


def score(
    left_primer: Oligo,
    right_primer: Oligo,
    amplicon: Span,
    amplicon_tm: float,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
    params: PrimerParameters,
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
        params: the set of primer3 parameters

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
        size_penalty = (amplicon.length - amplicon_sizes.opt) * params.amplicon_size_wt.gt
    else:
        size_penalty = (amplicon_sizes.opt - amplicon.length) * params.amplicon_size_wt.lt

    # The penalty for the amplicon melting temperature.
    # The difference in melting temperature between the calculated and optimal is weighted by the
    # product melting temperature.
    tm_penalty: float
    if amplicon_tms.opt == 0.0:
        tm_penalty = 0.0
    elif amplicon_tm > amplicon_tms.opt:
        tm_penalty = (amplicon_tm - amplicon_tms.opt) * params.amplicon_tm_wt.gt
    else:
        tm_penalty = (amplicon_tms.opt - amplicon_tm) * params.amplicon_tm_wt.lt

    # Put it all together
    return left_primer.penalty + right_primer.penalty + size_penalty + tm_penalty


def build_primer_pairs(  # noqa: C901
    left_primers: Sequence[Oligo],
    right_primers: Sequence[Oligo],
    target: Span,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
    max_heterodimer_tm: Optional[float],
    params: PrimerParameters,
    fasta_path: Path,
    thermo: Optional[Thermo] = None,
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
        params: the set of penalty params
        fasta_path: the path to the FASTA file from which the amplicon sequence will be retrieved.
        thermo: a [`Thermo`][prymer.Thermo] instance for performing thermodynamic calculations
            including amplicon tm; if not provided, a default Thermo instance will be created

    Returns:
        An iterator over all the valid primer pairs, sorted by primer pair penalty.
        Primer pairs with smaller penalties are returned first.
    """
    # Short circuit if we have no left primers or no right primers
    if not any(left_primers) or not any(right_primers):
        return

    if any(p.span.refname != target.refname for p in left_primers):
        raise ValueError("Left primers exist on different reference than target.")

    if any(p.span.refname != target.refname for p in right_primers):
        raise ValueError("Right primers exist on different reference than target.")

    # Sort the left and right primers
    left_primers = sorted(left_primers, key=lambda p: p.span.start)
    right_primers = sorted(right_primers, key=lambda p: p.span.end)

    # Grab the sequence we'll use to fill in the amplicon sequence
    with FastaFile(f"{fasta_path}") as fasta:
        region_start = min(p.span.start for p in left_primers)
        region_end = max(p.span.end for p in right_primers)
        bases = fasta.fetch(target.refname, region_start - 1, region_end)

    # Make sure we can do thermodynamic/Tm calculations
    thermo = thermo if thermo is not None else Thermo()

    # Each tuple is left_idx, right_idx, penalty, tm
    pairings: list[Tuple[int, int, float, float]] = []

    # generate all the primer pairs that don't violate hard size and Tm constraints
    first_right_primer_idx = 0

    # Nested loops over indices are used here so that we can skip potentially large chunks of
    # the cartesian product, based on the fact that we're sorting the left and right primers.
    # Two things are relied upon:
    #   1. If we encounter a left/right combo that either has the right primer leftward of the
    #      left primer _or_ generates a too-short amplicon, the neither that right primer nor
    #      any previous right primer can make a valid combination with any subsequent left primer.
    #   2. If we encounter a left/right combo that generates a too-large amplicon, then no
    #      subsequent right-primer can make a valid combination with that left primer
    for i in range(0, len(left_primers)):
        for j in range(first_right_primer_idx, len(right_primers)):
            lp = left_primers[i]
            rp = right_primers[j]

            # If the right primer isn't "to the right" of the left primer, move on
            if rp.span.start < lp.span.start or lp.span.end > rp.span.end:
                first_right_primer_idx = max(first_right_primer_idx, j + 1)
                continue

            amp_span = PrimerPair.calculate_amplicon_span(lp, rp)

            if amp_span.length < amplicon_sizes.min:
                first_right_primer_idx = max(first_right_primer_idx, j + 1)
                continue

            if amp_span.length > amplicon_sizes.max:
                break  # break in this case because all subsequent rps will yield longer amplicons

            # Since the amplicon span and the region_start are both 1-based, the minuend
            # becomes a zero-based offset
            amp_bases = bases[amp_span.start - region_start : amp_span.end - region_start + 1]
            amp_tm = thermo.tm(amp_bases)

            if amp_tm < amplicon_tms.min or amp_tm > amplicon_tms.max:
                continue

            penalty = score(
                left_primer=lp,
                right_primer=rp,
                amplicon=amp_span,
                amplicon_tm=amp_tm,
                amplicon_sizes=amplicon_sizes,
                amplicon_tms=amplicon_tms,
                params=params,
            )

            pairings.append((i, j, penalty, amp_tm))

    # Sort by the penalty, ascending
    pairings.sort(key=lambda tup: tup[2])

    for i, j, penalty, tm in pairings:
        lp = left_primers[i]
        rp = right_primers[j]

        if max_heterodimer_tm is not None:
            if thermo.heterodimer_tm(lp.bases, rp.bases) > max_heterodimer_tm:
                continue

        amp_bases = bases[lp.span.start - region_start : rp.span.end - region_start + 1]

        pp = PrimerPair(
            left_primer=lp,
            right_primer=rp,
            amplicon_sequence=amp_bases,
            amplicon_tm=tm,
            penalty=penalty,
        )

        yield pp
