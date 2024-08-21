"""
# Methods and class for building and scoring primer pairs from a list of left and right primers.

Typically, the first step is to build primer pairs from a list of left and right primers for a
given target with the [`build_primer_pairs()`][prymer.api.picking.build_primer_pairs]
method. The returned primer pairs are automatically scored (using the
[`score()`][prymer.api.picking.score] method), and returned sorted by the penalty
(increasing).

Next, the list of primer pairs for a given target are filtered based on various criteria using the
[`pick_top_primer_pairs()`][prymer.api.picking.pick_top_primer_pairs] method. Criteria
included but not limited to filtering for desired size and melting temperature ranges (see
[`is_acceptable_primer_pair()`][prymer.api.picking.is_acceptable_primer_pair]), not too
many off-targets, no self-dimers, and so on.  A given maximum number of passing primer pairs is
returned.

These two steps can be performed jointly using the
[`build_and_pick_primer_pairs()`][prymer.api.picking.build_and_pick_primer_pairs] method.

## Module contents

Contains the following public classes and methods:

- [`FilteringParams`][prymer.api.picking.FilteringParams] -- Stores the parameters used
    for filtering primers and primer pairs.
- [`score()`][prymer.api.picking.score] -- Scores the amplicon amplified by a primer pair
    in a manner similar to Primer3.
- [`check_primer_overlap()`][prymer.api.picking.check_primer_overlap] -- Checks that both
    (next) primers have at least `min_difference` bases that do not overlap the other (previous)
    primers.
- [`is_acceptable_primer_pair()`][prymer.api.picking.is_acceptable_primer_pair] -- Checks
    if a primer pair has the desired size range and melting temperature.
- [`build_primer_pairs()`][prymer.api.picking.build_primer_pairs] -- Builds primer pairs
    from individual left and primers.
- [`pick_top_primer_pairs()`][prymer.api.picking.pick_top_primer_pairs] -- Selects up to
    the given number of primer pairs from the given list of primer pairs.
- [`build_and_pick_primer_pairs()`][prymer.api.picking.build_and_pick_primer_pairs] --
    Builds primer pairs from individual left and primers and selects up to the given number of
    primer pairs from the given list of primer pairs.

"""

from dataclasses import dataclass
from typing import Callable
from typing import Iterable
from typing import Optional

from fgpyo.collections import PeekableIterator
from pysam import FastaFile

from prymer.api.melting import calculate_long_seq_tm
from prymer.api.minoptmax import MinOptMax
from prymer.api.primer import Primer
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.ntthal import NtThermoAlign
from prymer.offtarget.offtarget_detector import OffTargetDetector


@dataclass(frozen=True, init=True, slots=True)
class FilteringParams:
    """Stores the parameters used for filtering primers and primer pairs.

    Attributes:
        amplicon_sizes: the min, optimal, and max amplicon size
        amplicon_tms: the min, optimal, and max amplicon melting temperatures
        product_size_lt: penalty weight for products shorter than the optimal amplicon size
        product_size_gt: penalty weight for products longer than the optimal amplicon size
        product_tm_lt: penalty weight for products  with a Tm smaller the optimal amplicon Tm
        product_tm_gt: penalty weight for products  with a Tm larger the optimal amplicon Tm
        max_dimer_tm: the max primer dimer melting temperature allowed when filtering primer pairs
            for dimer formation.  Any primer pair with a dimer melting temperature above this
            threshold will be discarded.
        max_primer_hits: the maximum number of hits an individual primer can have in the genome
            before it is considered an invalid primer, and all primer pairs containing the primer
            failed.
        max_primer_pair_hits: the maximum number of amplicons a primer pair can make and be
            considered passing.
        three_prime_region_length: the number of bases at the 3' end of the primer in which the
            parameter `max_mismatches_in_three_prime_region` is evaluated.
        max_mismatches_in_three_prime_region: the maximum number of mismatches that are tolerated in
            the three prime region of each primer defined by `three_prime_region_length`.
        max_mismatches: the maximum number of mismatches that are tolerated in the full primer.
        min_primer_to_target_distance: minimum distance allowed between the end of a primer and
            start of the target -- aimed at centering the target in the amplicon.
        read_length: sequencing depth and maximum distance allowed between the start of a primer and
            the end of the target -- ensuring coverage during sequencing
        dist_from_primer_end_weight: weight determining importance of the distance between the
            primer and amplicon
        target_not_covered_by_read_length_weight: weight determining importance of coverage of the
            target by the sequencing depth
    """

    amplicon_sizes: MinOptMax[int]
    amplicon_tms: MinOptMax[float]
    product_size_lt: int
    product_size_gt: int
    product_tm_lt: float
    product_tm_gt: float
    max_dimer_tm: float = 55.0
    max_primer_hits: int = 500
    max_primer_pair_hits: int = 1
    three_prime_region_length: int = 1
    max_mismatches_in_three_prime_region: int = 2
    max_mismatches: int = 2
    min_primer_to_target_distance: int = 30
    read_length: int = 150
    dist_from_primer_end_weight: float = 100.0
    target_not_covered_by_read_length_weight: float = 100.0


def _dist_penalty(start: int, end: int, params: FilteringParams) -> float:
    """Returns a penalty that penalizes primers whose innermost base is closer than some
    minimum distance from the target.  The "distance" is zero if the innermost base overlaps the
    target, one if the innermost base abuts the target, etc.

    Args:
        start: the 1-based start position (inclusive)
        end: the 1-based end position (inclusive)
        params: the filtering parameters
    """
    diff = end - start
    if diff >= params.min_primer_to_target_distance:
        return 0.0
    else:
        return (params.min_primer_to_target_distance - diff) * params.dist_from_primer_end_weight


def _seq_penalty(start: int, end: int, params: FilteringParams) -> float:
    """Returns a penalty that penalizes amplicons where the target cannot be fully sequenced
    at the given read length starting from both ends of the amplicon.

    Args:
        start: the start coordinate of the span to consider
        end: the end coordinate of the span to consider
        params: the filtering parameters
    """

    amplicon_length = end - start + 1
    if amplicon_length <= params.read_length:
        return 0.0
    else:
        return (
            amplicon_length - params.read_length
        ) * params.target_not_covered_by_read_length_weight


def score(
    left: Primer,
    right: Primer,
    target: Span,
    amplicon: Span,
    amplicon_seq_or_tm: str | float,
    params: FilteringParams,
) -> float:
    """Score the amplicon in a manner similar to Primer3

    Sums the following:

    1. Primer penalties: the provided left and right primer penalties
    2. Amplicon size: The difference between the current amplicon size and optimal amplicon size
       scaled by the product size weight.  Is zero if the optimal amplicon size is zero.
    3. Amplicon Tm: The difference in melting temperature between the calculated and optimal,
       weighted by the product melting temperature.
    4. Inner primer distance: For primers whose innermost base is closer than some minimum distance
       from the target, the difference between the two distances scale by the corresponding weight.
    5. Amplicon length: The difference between the amplicon length and provided read length scaled
       by the corresponding weight.  Is zero when the amplicon is at most the read length.

    Args:
        left: the left primer
        right: the right primer
        target: the target mapping
        amplicon: the amplicon mapping
        amplicon_seq_or_tm: either the melting temperature of the amplicon, or the amplicon sequence
            from which the melting temperature of the amplicon will be calculated with
            `calculate_long_seq_tm`
        params: the filtering parameters

    Returns:
        the penalty for the whole amplicon.
    """
    # The penalty for the amplicon size:
    # 1. No penalty if the optimal amplicon size is zero
    # 2. The difference between the current amplicon size and optimal amplicon size scaled by the
    #    product size weight.  The product size weight is different depending on if the amplicon
    #    size is greater or less than the optimal amplicons.
    size_penalty: float
    if params.amplicon_sizes.opt == 0:
        size_penalty = 0.0
    elif amplicon.length > params.amplicon_sizes.opt:
        size_penalty = (amplicon.length - params.amplicon_sizes.opt) * params.product_size_gt
    else:
        size_penalty = (params.amplicon_sizes.opt - amplicon.length) * params.product_size_lt

    # The penalty for the amplicon melting temperature.
    # The difference in melting temperature between the calculated and optimal is weighted by the
    # product melting temperature.
    tm: float
    if isinstance(amplicon_seq_or_tm, str):
        tm = calculate_long_seq_tm(amplicon_seq_or_tm)
    else:
        tm = float(amplicon_seq_or_tm)
    tm_penalty: float
    if tm > params.amplicon_tms.opt:
        tm_penalty = (tm - params.amplicon_tms.opt) * params.product_tm_gt
    else:
        tm_penalty = (params.amplicon_tms.opt - tm) * params.product_tm_lt

    # Penalize primers whose innermost base is closer than some minimum distance from the target
    left_dist_penalty: float = _dist_penalty(start=left.span.end, end=target.start, params=params)
    right_dist_penalty: float = _dist_penalty(start=target.end, end=right.span.start, params=params)

    # Penalize amplicons where the target cannot be fully sequenced at the given read length
    # starting from both ends of the amplicon.
    left_seq_penalty: float = _seq_penalty(start=left.span.start, end=target.end, params=params)
    right_seq_penalty: float = _seq_penalty(start=target.start, end=right.span.end, params=params)

    # Put it all together
    return (
        left.penalty
        + right.penalty
        + size_penalty
        + tm_penalty
        + left_dist_penalty
        + right_dist_penalty
        + left_seq_penalty
        + right_seq_penalty
    )


def check_primer_overlap(
    prev_left: Span,
    prev_right: Span,
    next_left: Span,
    next_right: Span,
    min_difference: int,
) -> bool:
    """Check that both (next) primers have at least `min_difference` bases that do not overlap the
     other (previous) primers.

    Args:
        prev_left: left primer mapping from the previous primer pair.
        prev_right: right primer mapping from the previous primer pair
        next_left: left primer mapping from the current primer pair being compared
        next_right: right primer mapping from the current primer pair being compared
        min_difference: the minimum number of bases that must differ between two primers

    Returns:
        True if the primers are sufficiently different
    """
    left = prev_left.length_of_overlap_with(next_left) + min_difference
    right = prev_right.length_of_overlap_with(next_right) + min_difference
    return (
        prev_left.length >= left
        and next_left.length >= left
        and prev_right.length >= right
        and next_right.length >= right
    )


def is_acceptable_primer_pair(primer_pair: PrimerPair, params: FilteringParams) -> bool:
    """Determine if a primer pair can be kept considering:

    1. the primer pair is in the desired size range
    2. The primer pair has a melting temperature in the desired range

    Args:
        primer_pair: the primer pair.
        params: the parameters used for filtering

    Returns:
       True if the primer pair passes initial checks and can be considered later on.
    """
    return (
        params.amplicon_sizes.min <= primer_pair.amplicon.length <= params.amplicon_sizes.max
        and params.amplicon_tms.min <= primer_pair.amplicon_tm <= params.amplicon_tms.max
    )


def build_primer_pairs(
    lefts: Iterable[Primer],
    rights: Iterable[Primer],
    target: Span,
    params: FilteringParams,
    fasta: FastaFile,
) -> list[PrimerPair]:
    """Builds primer pairs from individual left and primers.

    Args:
        lefts: the left primers
        rights: the right primers
        target: the genome mapping for the target
        params: the parameters used for filtering
        fasta: the FASTA file from which the amplicon sequence will be retrieved.

    Returns:
        the list of primer pairs, sorted by penalty (increasing)
    """
    # generate all the primer pairs
    primer_pairs: list[PrimerPair] = []
    for left in lefts:
        for right in rights:
            if left.span.refname != right.span.refname:
                raise ValueError(
                    "Cannot create a primer pair from left and right primers on different"
                    f"references: left: '{left.span.refname}' right: '{right.span.refname}'"
                )
            amplicon_mapping = Span(
                refname=target.refname, start=left.span.start, end=right.span.end
            )
            amplicon_bed = amplicon_mapping.get_bedlike_coords()  # since fasta.fetch is 0-based
            amplicon_sequence = fasta.fetch(
                reference=target.refname, start=amplicon_bed.start, end=amplicon_bed.end
            )
            amplicon_penalty = score(
                left=left,
                right=right,
                target=target,
                amplicon=amplicon_mapping,
                amplicon_seq_or_tm=amplicon_sequence,
                params=params,
            )
            pp = PrimerPair(
                left_primer=left,
                right_primer=right,
                amplicon_sequence=amplicon_sequence,
                amplicon_tm=calculate_long_seq_tm(amplicon_sequence),
                penalty=amplicon_penalty,
            )
            primer_pairs.append(pp)

    # sort by smallest penalty
    return sorted(primer_pairs, key=lambda pp: pp.penalty)


def pick_top_primer_pairs(
    primer_pairs: list[PrimerPair],
    num_primers: int,
    min_difference: int,
    params: FilteringParams,
    offtarget_detector: OffTargetDetector,
    is_dimer_tm_ok: Callable[[str, str], bool],
) -> list[PrimerPair]:
    """Selects up to the given number of primer pairs from the given list of primer pairs.

    The primer pairs are selected in the order of lowest penalty.

    The primer pairs must:
    1. Have an amplicon in desired size range (see `is_acceptable_primer_pair`).
    2. Have an amplicon melting temperature in the desired range (see `is_acceptable_primer_pair`).
    3. Not have too many off-targets (see `OffTargetDetector#check_one()`).
    4. Not have primer pairs that overlap too much (see `check_primer_overlap`).
    5. Not form a dimer with a melting temperature above a specified threshold (see `max_dimer_tm`).

    The _order_ of these checks are important as some of them are expensive to compute.

    Args:
        primer_pairs: the list of all primer pairs.
        num_primers: the number of primer pairs to return for the target
        min_difference: the minimum base difference between two primers that we will tolerate.
        params: the parameters used for filtering
        offtarget_detector: the off-target detector.
        is_dimer_tm_ok: a function to use for checking for dimers. Should return true if the pair
            passes (e.g. the dimer check)

    Returns:
        Up to `num_primers` primer pairs
    """
    selected: list[PrimerPair] = []
    pp_iter = PeekableIterator(primer_pairs)
    last_pp: Optional[PrimerPair] = None
    while len(selected) < num_primers and pp_iter.can_peek():
        pp = next(pp_iter)

        # Enforce that primers were ordered by penalty (larger value first)
        if last_pp is not None and pp.penalty < last_pp.penalty:
            raise ValueError("Primers must be given in order by penalty (increasing value).")
        last_pp = pp

        # Check some basic properties
        if not is_acceptable_primer_pair(primer_pair=pp, params=params):
            continue

        # Check for off-targets
        if not offtarget_detector.check_one(primer_pair=pp).passes:
            continue

        # Check that both primers have at least `min_difference` bases that do not overlap the
        # other primer.
        okay = all(
            check_primer_overlap(
                prev.left_primer.span,
                prev.right_primer.span,
                pp.left_primer.span,
                pp.right_primer.span,
                min_difference=min_difference,
            )
            for prev in selected
        )
        if not okay:
            continue

        # Check dimer Tm
        if not is_dimer_tm_ok(pp.left_primer.bases, pp.right_primer.bases):
            continue

        selected.append(pp)

    return selected


def build_and_pick_primer_pairs(
    lefts: Iterable[Primer],
    rights: Iterable[Primer],
    target: Span,
    num_primers: int,
    min_difference: int,
    params: FilteringParams,
    offtarget_detector: OffTargetDetector,
    dimer_checker: NtThermoAlign,
    fasta: FastaFile,
) -> list[PrimerPair]:
    """Picks up to `num_primers` primer pairs.

    Args:
        lefts: the left primers
        rights: the right primers
        target: the genome mapping for the target
        num_primers: the number of primer pairs to return for the target.
        min_difference: the minimum base difference between two primers that we will tolerate.
        params: the parameters used for filtering.
        offtarget_detector: the off-target detector.
        dimer_checker: the primer-dimer melting temperature checker.
        fasta: the FASTA file from which the amplicon sequence will be retrieved.

    Returns:
        the list of primer pairs, sorted by penalty (increasing)
    """
    # build the list of primer pairs
    primer_pairs = build_primer_pairs(
        lefts=lefts, rights=rights, target=target, params=params, fasta=fasta
    )

    # select the primer pairs
    selected: list[PrimerPair] = pick_top_primer_pairs(
        primer_pairs=primer_pairs,
        num_primers=num_primers,
        min_difference=min_difference,
        params=params,
        offtarget_detector=offtarget_detector,
        is_dimer_tm_ok=lambda s1, s2: (
            dimer_checker.duplex_tm(s1=s1, s2=s2) <= params.max_dimer_tm
        ),
    )

    return selected
