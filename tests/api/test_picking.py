from collections import Counter
from dataclasses import dataclass
from dataclasses import replace
from pathlib import Path
from typing import Any
from typing import Optional
from typing import Tuple

import pysam
import pytest
from fgpyo.sequence import reverse_complement

from prymer.api import FilteringParams
from prymer.api import MinOptMax
from prymer.api import Oligo
from prymer.api import PrimerPair
from prymer.api import Span
from prymer.api import Strand
from prymer.api import build_and_pick_primer_pairs
from prymer.api import build_primer_pairs
from prymer.api import pick_top_primer_pairs
from prymer.api.melting import calculate_long_seq_tm
from prymer.api.picking import _dist_penalty
from prymer.api.picking import _seq_penalty
from prymer.api.picking import check_primer_overlap
from prymer.api.picking import is_acceptable_primer_pair
from prymer.api.picking import score as picking_score
from prymer.ntthal import NtThermoAlign
from prymer.offtarget import OffTargetDetector


@pytest.fixture
def filter_params() -> FilteringParams:
    return FilteringParams(
        amplicon_sizes=MinOptMax(100, 150, 200),
        amplicon_tms=MinOptMax(50.0, 55.0, 60.0),
        product_size_lt=1,
        product_size_gt=1,
        product_tm_lt=0.0,
        product_tm_gt=0.0,
        min_primer_to_target_distance=0,
        read_length=100000,
    )


@pytest.mark.parametrize(
    "start, end, min_primer_to_target_distance, dist_from_primer_end_weight, penalty",
    [
        (10, 20, 10, 2.0, 0.0),  # 10bp away, at min distance
        (10, 19, 10, 2.0, 2.0),  # 9bp away, below min distance
        (10, 21, 10, 2.0, 0.0),  # 11bp away, beyond min distance
        (10, 22, 32, 3.0, 60.0),  # 2bp away, one over min distance, larger weight
        (10, 22, 32, 5.0, 100.0),  # 2bp away, one over min distance, large weight
    ],
)
def test_dist_penalty(
    start: int,
    end: int,
    min_primer_to_target_distance: int,
    dist_from_primer_end_weight: float,
    penalty: float,
    filter_params: FilteringParams,
) -> None:
    params = replace(
        filter_params,
        min_primer_to_target_distance=min_primer_to_target_distance,
        dist_from_primer_end_weight=dist_from_primer_end_weight,
    )
    assert _dist_penalty(start=start, end=end, params=params) == penalty


@pytest.mark.parametrize(
    "start, end, read_length, target_not_covered_by_read_length_weight, penalty",
    [
        (10, 10, 1, 2.0, 0.0),  # 1bp amplicon_length length, at read length
        (10, 10, 2, 2.0, 0.0),  # 1bp amplicon_length length, one shorter than read length
        (10, 10, 0, 2.0, 2.0),  # 1bp amplicon_length length, one longer than read length
        (10, 100, 91, 2.0, 0.0),  # 91bp amplicon length, at read length
        (10, 100, 90, 2.0, 2.0),  # 1bp amplicon_length length, one longer than read length
        (10, 100, 91, 4.0, 0.0),  # 1bp amplicon_length length, one below minimum, larger weight
        (10, 100, 50, 5.0, 205.0),  # 1bp amplicon_length length, far below minimum, large weight
    ],
)
def test_seq_penalty(
    start: int,
    end: int,
    read_length: int,
    target_not_covered_by_read_length_weight: float,
    penalty: float,
    filter_params: FilteringParams,
) -> None:
    params = replace(
        filter_params,
        read_length=read_length,
        target_not_covered_by_read_length_weight=target_not_covered_by_read_length_weight,
    )
    assert _seq_penalty(start=start, end=end, params=params) == penalty


def build_primer_pair(amplicon_length: int, tm: float) -> PrimerPair:
    left_primer = Oligo(
        tm=0,
        penalty=0,
        span=Span(refname="1", start=1, end=max(1, amplicon_length // 4)),
    )
    right_primer = Oligo(
        tm=0,
        penalty=0,
        span=Span(
            refname="1", start=max(1, amplicon_length - (amplicon_length // 4)), end=amplicon_length
        ),
    )
    return PrimerPair(
        left_primer=left_primer,
        right_primer=right_primer,
        amplicon_tm=tm,
        penalty=0,
    )


@pytest.mark.parametrize(
    "prev_left, prev_right, next_left, next_right, min_difference, expected",
    [
        # no overlap -> True
        (
            Span("chr1", 100, 120),
            Span("chr1", 200, 220),
            Span("chr1", 120, 140),
            Span("chr1", 220, 240),
            10,
            True,
        ),
        # left at min_difference -> True
        (
            Span("chr1", 100, 120),
            Span("chr1", 200, 220),
            Span("chr1", 110, 130),
            Span("chr1", 220, 240),
            10,
            True,
        ),
        # left one less than min_difference -> False
        (
            Span("chr1", 100, 120),
            Span("chr1", 200, 220),
            Span("chr1", 109, 129),
            Span("chr1", 220, 240),
            10,
            False,
        ),
        # right at min_difference -> True
        (
            Span("chr1", 100, 120),
            Span("chr1", 200, 220),
            Span("chr1", 120, 140),
            Span("chr1", 210, 230),
            10,
            True,
        ),
        # right one less than min_difference -> False
        (
            Span("chr1", 100, 120),
            Span("chr1", 200, 220),
            Span("chr1", 120, 140),
            Span("chr1", 209, 229),
            10,
            False,
        ),
        # both at min_difference -> True
        (
            Span("chr1", 100, 120),
            Span("chr1", 200, 220),
            Span("chr1", 110, 130),
            Span("chr1", 210, 230),
            10,
            True,
        ),
        # both one less than min_difference -> False
        (
            Span("chr1", 100, 120),
            Span("chr1", 200, 220),
            Span("chr1", 109, 129),
            Span("chr1", 209, 229),
            10,
            False,
        ),
    ],
)
def test_check_primer_overlap(
    prev_left: Span,
    prev_right: Span,
    next_left: Span,
    next_right: Span,
    min_difference: int,
    expected: bool,
) -> None:
    assert (
        check_primer_overlap(
            prev_left=prev_left,
            prev_right=prev_right,
            next_left=next_left,
            next_right=next_right,
            min_difference=min_difference,
        )
        == expected
    )
    assert (
        check_primer_overlap(
            prev_left=next_left,
            prev_right=next_right,
            next_left=prev_left,
            next_right=prev_right,
            min_difference=min_difference,
        )
        == expected
    )


@pytest.mark.parametrize(
    "pair, expected",
    [
        (build_primer_pair(amplicon_length=150, tm=55), True),  # OK at optimal
        (build_primer_pair(amplicon_length=100, tm=55), True),  # OK, amplicon at lower boundary
        (build_primer_pair(amplicon_length=99, tm=55), False),  # NOK, amplicon < lower boundary
        (build_primer_pair(amplicon_length=200, tm=55), True),  # OK, amplicon at upper boundary
        (build_primer_pair(amplicon_length=201, tm=55), False),  # NOK, amplicon > upper boundary
        (build_primer_pair(amplicon_length=150, tm=50), True),  # OK, tm at lower boundary
        (build_primer_pair(amplicon_length=150, tm=49), False),  # NOK, tm < lower boundary
        (build_primer_pair(amplicon_length=150, tm=60), True),  # OK, tm at upper boundary
        (build_primer_pair(amplicon_length=150, tm=61), False),  # NOK, tm > upper boundary
    ],
)
def test_is_acceptable_primer_pair(pair: PrimerPair, expected: bool) -> None:
    params = FilteringParams(
        amplicon_sizes=MinOptMax(100, 150, 200),
        amplicon_tms=MinOptMax(50.0, 55.0, 60.0),
        product_size_lt=1,
        product_size_gt=1,
        product_tm_lt=0.0,
        product_tm_gt=0.0,
        min_primer_to_target_distance=0,
        read_length=100000,
    )
    assert is_acceptable_primer_pair(primer_pair=pair, params=params) == expected


@dataclass(init=True, frozen=True)
class ScoreInput:
    left: Oligo
    right: Oligo
    target: Span
    amplicon: Span
    amplicon_sequence: str

    def __post_init__(self) -> None:
        primer_pair = PrimerPair(
            left_primer=self.left,
            right_primer=self.right,
            amplicon_tm=0,
            penalty=0,
        )
        object.__setattr__(self, "primer_pair", primer_pair)


def _score_input() -> ScoreInput:
    l_mapping = Span(refname="1", start=95, end=125)
    r_mapping = Span(refname="1", start=195, end=225)
    amplicon = Span(refname="1", start=l_mapping.end, end=r_mapping.start)
    target = Span(refname="1", start=l_mapping.end + 10, end=r_mapping.start - 20)
    return ScoreInput(
        left=Oligo(penalty=0, tm=0, span=l_mapping),
        right=Oligo(penalty=0, tm=0, span=r_mapping),
        target=target,
        amplicon=amplicon,
        amplicon_sequence="A" * amplicon.length,
    )


@pytest.fixture()
def score_input() -> ScoreInput:
    return _score_input()


# Developer Note: for `score()`, we must have:
# 1. params.amplicon_sizes.opt == 0
# 2. params.amplicon_tms.opt is exactly the _calculate_long_seq_tm of the amplicon sequence
# 3. filtering_params.min_primer_to_target_distance is zero,
#    (target.start >= left.span.end)
# 4. filtering_params.min_primer_to_target_distance is zero,
#    (right.span.start >= target.end)
# 5. filtering_params.read_length is large, (target.end >= left.span.start)
# 6. filtering_params.read_length is large, (right.span.end >= target.start)
def _zero_score_filtering_params(score_input: ScoreInput) -> FilteringParams:
    amplicon_tm = calculate_long_seq_tm(score_input.amplicon_sequence)
    return FilteringParams(
        # for size_penalty to be zero
        amplicon_sizes=MinOptMax(0, 0, 200),
        product_size_gt=0,
        product_size_lt=0,
        # for tm_penalty to be zero
        amplicon_tms=MinOptMax(amplicon_tm, amplicon_tm, amplicon_tm),
        product_tm_gt=0,
        product_tm_lt=0,
        # for left_dist_penalty and left_dist_penalty to be zero
        min_primer_to_target_distance=0,
        dist_from_primer_end_weight=0,
        # for left_seq_penalty and left_seq_penalty to be zero
        read_length=100000,
        target_not_covered_by_read_length_weight=0,
    )


@pytest.fixture()
def zero_score_filtering_params(score_input: ScoreInput) -> FilteringParams:
    return _zero_score_filtering_params(score_input=score_input)


def test_zero_score(
    score_input: ScoreInput,
    zero_score_filtering_params: FilteringParams,
) -> None:
    assert (
        picking_score(
            left=score_input.left,
            right=score_input.right,
            target=score_input.target,
            amplicon=score_input.amplicon,
            amplicon_seq_or_tm=score_input.amplicon_sequence,
            params=zero_score_filtering_params,
        )
        == 0.0
    )


def test_zero_score_with_amplicon_tm(
    score_input: ScoreInput,
    zero_score_filtering_params: FilteringParams,
) -> None:
    amplicon_tm: float = calculate_long_seq_tm(score_input.amplicon_sequence)
    assert (
        picking_score(
            left=score_input.left,
            right=score_input.right,
            target=score_input.target,
            amplicon=score_input.amplicon,
            amplicon_seq_or_tm=amplicon_tm,
            params=zero_score_filtering_params,
        )
        == 0.0
    )


# Notes on score_input.
# - amplicon length is 71bp (for size_penalty)
# - amplicon Tm is 66.30319122042972 (for tm_penalty)
# - primer to target distance is 10 for left_dist_penalty and 20 right_dist_penalty.  This also
#   means we need a read length of at least 81bp (=71 + 10) for the left_seq_penalty to be zero
#   and at least 91bp (=71 + 20) for the right_seq_penalty to be zero
@pytest.mark.parametrize(
    "kwargs, expected_score",
    [
        # size_penalty: amplicon size is too large (71bp > 61bp, diff = 10, so score = (71-61)*5)
        ({"product_size_gt": 5.0, "amplicon_sizes": MinOptMax(61, 61, 61)}, (71 - 61) * 5),
        # size_penalty: amplicon size is too small (71bp < 81bp, diff = 10, so score = (71-61)*4)
        ({"product_size_lt": 4.0, "amplicon_sizes": MinOptMax(81, 81, 81)}, (71 - 61) * 4),
        # tm_penalty: amplicon Tm is too large
        (
            {"product_tm_gt": 5.0, "amplicon_tms": MinOptMax(60, 60, 60)},
            5.0 * (66.30319122042972 - 60),
        ),
        # tm_penalty: amplicon Tm is too small
        (
            {"product_tm_lt": 4.0, "amplicon_tms": MinOptMax(70, 70, 70)},
            4.0 * (70 - 66.30319122042972),
        ),
        # both [left/right]_dist_penalty are zero: diff > min_primer_to_target_distance
        ({"min_primer_to_target_distance": 10, "dist_from_primer_end_weight": 5.0}, 0),
        # left_dist_penalty: diff < min_primer_to_target_distance
        ({"min_primer_to_target_distance": 11, "dist_from_primer_end_weight": 5.0}, (1 * 5.0)),
        ({"min_primer_to_target_distance": 20, "dist_from_primer_end_weight": 5.0}, (10 * 5.0)),
        # both [left/right]_dist_penalty: diff < min_primer_to_target_distance, 11bp for left, and
        # 1bp for right, so 12bp overall
        ({"min_primer_to_target_distance": 21, "dist_from_primer_end_weight": 5.0}, (12 * 5.0)),
        # both[left/right]_seq_penalty: zero since read length >= 81bp (left) and >= 91bp (right)
        ({"read_length": 91, "target_not_covered_by_read_length_weight": 10.0}, 0),
        # right_seq_penalty: less than 91bp for the right but at 81bp for the left
        ({"read_length": 81, "target_not_covered_by_read_length_weight": 10.0}, (10 * 10)),
        # both[left/right]_seq_penalty: zero since read length <>=> 81bp (left) and < 91bp (right)
        ({"read_length": 71, "target_not_covered_by_read_length_weight": 10.0}, 10.0 * (10 + 20)),
    ],
)
def test_score(
    score_input: ScoreInput,
    zero_score_filtering_params: FilteringParams,
    kwargs: dict[str, Any],
    expected_score: float,
) -> None:
    params = replace(zero_score_filtering_params, **kwargs)
    assert (
        picking_score(
            left=score_input.left,
            right=score_input.right,
            target=score_input.target,
            amplicon=score_input.amplicon,
            amplicon_seq_or_tm=score_input.amplicon_sequence,
            params=params,
        )
        == expected_score
    )


@pytest.mark.parametrize(
    "left_penalty, right_penalty", [(0, 0), (10.0, 0.0), (0, 12.0), (13.0, 14.0)]
)
def test_score_primer_primer_penalties(
    score_input: ScoreInput,
    zero_score_filtering_params: FilteringParams,
    left_penalty: float,
    right_penalty: float,
) -> None:
    left = replace(score_input.left, penalty=left_penalty)
    right = replace(score_input.right, penalty=right_penalty)
    assert picking_score(
        left=left,
        right=right,
        target=score_input.target,
        amplicon=score_input.amplicon,
        amplicon_seq_or_tm=score_input.amplicon_sequence,
        params=zero_score_filtering_params,
    ) == (left_penalty + right_penalty)


def test_primer_pairs(
    score_input: ScoreInput, zero_score_filtering_params: FilteringParams, genome_ref: Path
) -> None:
    primer_length: int = 30
    target = Span(refname="chr1", start=100, end=250)

    # tile some left primers
    lefts = []
    rights = []
    for offset in range(0, 50, 5):
        # left
        left_end = target.start - offset
        left_start = left_end - primer_length
        left = Oligo(
            penalty=-offset, tm=0, span=Span(refname=target.refname, start=left_start, end=left_end)
        )
        lefts.append(left)
        # right
        right_start = target.end + offset
        right_end = right_start + primer_length
        right = Oligo(
            penalty=offset,
            tm=0,
            span=Span(refname=target.refname, start=right_start, end=right_end),
        )
        rights.append(right)

    with pysam.FastaFile(f"{genome_ref}") as fasta:
        primer_pairs = build_primer_pairs(
            lefts=lefts,
            rights=rights,
            target=target,
            params=zero_score_filtering_params,
            fasta=fasta,
        )
        assert len(primer_pairs) == len(lefts) * len(rights)
        last_penalty = primer_pairs[0].penalty
        primer_counter: Counter[Oligo] = Counter()
        for pp in primer_pairs:
            assert pp.left_primer in lefts
            assert pp.right_primer in rights
            primer_counter[pp.left_primer] += 1
            primer_counter[pp.right_primer] += 1
            # by design, only the left/right penalties contribute to the primer pair penalty
            assert pp.penalty == pp.left_primer.penalty + pp.right_primer.penalty
            # at least check that the amplicon Tm is non-zero
            assert pp.amplicon_tm > 0
            # at least check that the amplicon sequence retrieved is the correct length
            assert len(pp.amplicon_sequence) == pp.amplicon.length
            # check that the primer pairs are sorted by increasing penalty
            assert last_penalty <= pp.penalty
            last_penalty = pp.penalty
        # make sure we see all the primers the same # of times!
        items = primer_counter.items()
        assert len(set(i[0] for i in items)) == len(lefts) + len(rights)  # same primers
        assert len(set(i[1] for i in items)) == 1  # same counts for each primer


def test_primer_pairs_except_different_references(
    score_input: ScoreInput, zero_score_filtering_params: FilteringParams, genome_ref: Path
) -> None:
    # all primers (both left and right) _should_ be on the same reference (refname).  Add a test
    # that raises an exception (in PrimerPair) that doesn't.
    with pysam.FastaFile(f"{genome_ref}") as fasta:
        with pytest.raises(ValueError, match="Cannot create a primer pair"):
            # change the reference for the right primer
            right = replace(score_input.right, span=Span(refname="Y", start=195, end=225))
            build_primer_pairs(
                lefts=[score_input.left],
                rights=[right],
                target=score_input.target,
                params=zero_score_filtering_params,
                fasta=fasta,
            )


def _picking_ref() -> Path:
    return Path(__file__).parent / "data" / "picking.fa"


@pytest.fixture(scope="session")
def picking_ref() -> Path:
    return _picking_ref()


def _target() -> Span:
    target: Span = Span(refname="chr1", start=150, end=250)
    return target


def _primer_pair(
    target_offset: int = 0,
    primer_length: int = 30,
    amplicon_tm: Optional[float] = None,
    penalty: int = 0,
    target: Optional[Span] = None,
    params: Optional[FilteringParams] = None,
) -> PrimerPair:
    if target is None:
        target = _target()
    if params is None:
        params = _zero_score_filtering_params(score_input=_score_input())
    fasta = pysam.FastaFile(f"{_picking_ref()}")
    if amplicon_tm is None:
        amplicon_tm = params.amplicon_tms.opt
    left_span = Span(
        refname=target.refname,
        start=target.start - primer_length - target_offset + 1,  # +1 for 1-based inclusive
        end=target.start - target_offset,
        strand=Strand.POSITIVE,
    )
    assert left_span.length == primer_length
    left_bases = fasta.fetch(
        reference=left_span.refname, start=left_span.start - 1, end=left_span.end
    )
    right_span = Span(
        refname=target.refname,
        start=target.end + target_offset,
        end=target.end + primer_length + target_offset - 1,  # -1 for 1-based inclusive
        strand=Strand.NEGATIVE,
    )
    assert right_span.length == primer_length
    right_bases = fasta.fetch(
        reference=left_span.refname, start=right_span.start - 1, end=right_span.end
    )
    right_bases = reverse_complement(right_bases)
    fasta.close()
    return PrimerPair(
        left_primer=Oligo(bases=left_bases, penalty=0, tm=0, span=left_span),
        right_primer=Oligo(bases=right_bases, penalty=0, tm=0, span=right_span),
        amplicon_tm=amplicon_tm,
        penalty=penalty,
    )


def _pick_top_primer_pairs(
    params: FilteringParams,
    picking_ref: Path,
    primer_pairs: list[PrimerPair],
    max_primer_hits: int,
    max_primer_pair_hits: int,
    min_difference: int = 1,
) -> list[PrimerPair]:
    offtarget_detector = OffTargetDetector(
        ref=picking_ref,
        max_primer_hits=max_primer_hits,
        max_primer_pair_hits=max_primer_pair_hits,
        three_prime_region_length=5,
        max_mismatches_in_three_prime_region=0,
        max_mismatches=0,
        max_amplicon_size=params.amplicon_sizes.max,
    )
    dimer_checker = NtThermoAlign()

    picked = pick_top_primer_pairs(
        primer_pairs=primer_pairs,
        num_primers=len(primer_pairs),
        min_difference=min_difference,
        params=params,
        offtarget_detector=offtarget_detector,
        is_dimer_tm_ok=lambda s1, s2: (
            dimer_checker.duplex_tm(s1=s1, s2=s2) <= params.max_dimer_tm
        ),
    )
    offtarget_detector.close()
    dimer_checker.close()

    return picked


_PARAMS: FilteringParams = _zero_score_filtering_params(_score_input())
"""Filter parameters for use in creating the test cases for `test_pick_top_primer_pairs`"""


@pytest.mark.parametrize(
    "picked, primer_pair",
    [
        # primer pair passes all primer/primer-pair-specific filters
        (True, _primer_pair()),
        # too small amplicon size
        (False, _primer_pair(amplicon_tm=_PARAMS.amplicon_sizes.min - 1)),
        # too big amplicon size
        (False, _primer_pair(amplicon_tm=_PARAMS.amplicon_sizes.max + 1)),
        # too low amplicon Tm
        (False, _primer_pair(amplicon_tm=_PARAMS.amplicon_tms.min - 1)),
        # too large amplicon Tm
        (False, _primer_pair(amplicon_tm=_PARAMS.amplicon_tms.max + 1)),
    ],
)
def test_pick_top_primer_pairs_individual_primers(
    picking_ref: Path,
    zero_score_filtering_params: FilteringParams,
    picked: bool,
    primer_pair: PrimerPair,
) -> None:
    expected = [primer_pair] if picked else []
    assert expected == _pick_top_primer_pairs(
        params=zero_score_filtering_params,
        picking_ref=picking_ref,
        primer_pairs=[primer_pair],
        max_primer_hits=2,
        max_primer_pair_hits=2,
    )


def test_pick_top_primer_pairs_dimer_tm_too_large(
    picking_ref: Path,
    zero_score_filtering_params: FilteringParams,
) -> None:
    pp = _primer_pair()

    # get the Tm for the primer pair
    duplex_tm = NtThermoAlign().duplex_tm(s1=pp.left_primer.bases, s2=pp.right_primer.bases)

    # at the maximum dimer Tm
    params = replace(zero_score_filtering_params, max_dimer_tm=duplex_tm)
    assert [pp] == _pick_top_primer_pairs(
        params=params,
        picking_ref=picking_ref,
        primer_pairs=[pp],
        max_primer_hits=2,
        max_primer_pair_hits=2,
    )

    # **just** over the maximum dimer Tm
    params = replace(zero_score_filtering_params, max_dimer_tm=duplex_tm - 0.0001)
    assert [] == _pick_top_primer_pairs(
        params=params,
        picking_ref=picking_ref,
        primer_pairs=[pp],
        max_primer_hits=2,
        max_primer_pair_hits=2,
    )


_PRIMER_PAIRS_PICKING: list[PrimerPair] = [
    _primer_pair(target_offset=10, penalty=2),  # left.start=111
    _primer_pair(target_offset=6, penalty=2),  # 4bp from the previous, left.start=115
    _primer_pair(target_offset=3, penalty=3),  # 3bp from the previous, left.start=118
    _primer_pair(target_offset=1, penalty=4),  # 2bp from the previous, left.start=120
    _primer_pair(target_offset=0, penalty=5),  # 1bp from the previous, left.start=121
    _primer_pair(target_offset=0, penalty=6),  # same as the previous, left.start=121
]


@dataclass(frozen=True)
class _PickingTestCase:
    min_difference: int
    primer_pairs: list[Tuple[bool, PrimerPair]]


@pytest.mark.parametrize(
    "test_case",
    [
        # primer pairs that overlap each other fully, so second one is discarded
        _PickingTestCase(
            min_difference=1,
            primer_pairs=[
                (True, _primer_pair(target_offset=1, penalty=2)),  # first primer pair, kept
                (
                    False,
                    _primer_pair(target_offset=1, penalty=2),
                ),  # same as the previous, discarded
            ],
        ),
        # primer pairs that are offset by 1bp (equal to min-difference) so both are kept
        _PickingTestCase(
            min_difference=1,
            primer_pairs=[
                (True, _primer_pair(target_offset=0, penalty=2)),
                (True, _primer_pair(target_offset=1, penalty=2)),
            ],
        ),
        # primer pairs that are offset by 1bp (less than min-difference) so the first is kept
        _PickingTestCase(
            min_difference=2,
            primer_pairs=[
                (True, _primer_pair(target_offset=0, penalty=2)),
                (False, _primer_pair(target_offset=1, penalty=2)),
            ],
        ),
        _PickingTestCase(
            min_difference=3,
            primer_pairs=list(
                zip(
                    # 111,  115,  118,   120,  121,   121
                    [True, True, True, False, True, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
        _PickingTestCase(
            min_difference=4,
            primer_pairs=list(
                zip(
                    # 111,  115,   118,  120,   121,   121
                    [True, True, False, True, False, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
        _PickingTestCase(
            min_difference=5,
            primer_pairs=list(
                zip(
                    # 111,   115,  118,   120,   121,   121
                    [True, False, True, False, False, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
        _PickingTestCase(
            min_difference=6,
            primer_pairs=list(
                zip(
                    # 111,   115,  118,   120,   121,   121
                    [True, False, True, False, False, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
        _PickingTestCase(
            min_difference=7,
            primer_pairs=list(
                zip(
                    # 111,   115,  118,   120,   121,   121
                    [True, False, True, False, False, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
        _PickingTestCase(
            min_difference=8,
            primer_pairs=list(
                zip(
                    # 111,   115,   118,  120,   121,   121
                    [True, False, False, True, False, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
        _PickingTestCase(
            min_difference=9,
            primer_pairs=list(
                zip(
                    # 111,   115,   118,  120,   121,   121
                    [True, False, False, True, False, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
        _PickingTestCase(
            min_difference=10,
            primer_pairs=list(
                zip(
                    # 111,   115,   118,   120,  121,   121
                    [True, False, False, False, True, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
        _PickingTestCase(
            min_difference=11,
            primer_pairs=list(
                zip(
                    # 111,   115,   118,   120,   121,   121
                    [True, False, False, False, False, False],
                    _PRIMER_PAIRS_PICKING,
                    strict=True,
                )
            ),
        ),
    ],
)
def test_pick_top_primer_pairs_multiple_primers(
    picking_ref: Path,
    zero_score_filtering_params: FilteringParams,
    test_case: _PickingTestCase,
) -> None:
    in_primer_pairs = [pp[1] for pp in test_case.primer_pairs]

    # First check that picked _alone_ they will be kept
    for primer_pair in in_primer_pairs:
        assert [primer_pair] == _pick_top_primer_pairs(
            params=zero_score_filtering_params,
            picking_ref=picking_ref,
            primer_pairs=[primer_pair],
            max_primer_hits=2,
            max_primer_pair_hits=2,
            min_difference=test_case.min_difference,
        )

    # Next, check that picked _together_ only those primer pairs with "True" associated with them
    # are kept
    expected = [pp[1] for pp in test_case.primer_pairs if pp[0]]
    assert expected == _pick_top_primer_pairs(
        params=zero_score_filtering_params,
        picking_ref=picking_ref,
        primer_pairs=in_primer_pairs,
        max_primer_hits=2,
        max_primer_pair_hits=2,
        min_difference=test_case.min_difference,
    )


def test_pick_top_primer_pairs_input_order(
    picking_ref: Path,
    zero_score_filtering_params: FilteringParams,
) -> None:
    with pytest.raises(ValueError, match="Primers must be given in order by penalty"):
        _pick_top_primer_pairs(
            params=zero_score_filtering_params,
            picking_ref=picking_ref,
            primer_pairs=[_primer_pair(penalty=6), _primer_pair(penalty=5)],
            max_primer_hits=2,
            max_primer_pair_hits=2,
            min_difference=0,
        )


def test_pick_top_primer_pairs_too_many_off_targets(
    picking_ref: Path,
    zero_score_filtering_params: FilteringParams,
) -> None:
    """The contig is duplicated wholesale in `picking.fa`, so we should get two hits per primer.
    Thus, setting `max_primer_hits=1` will return "too many hits"""
    assert [] == _pick_top_primer_pairs(
        params=zero_score_filtering_params,
        picking_ref=picking_ref,
        primer_pairs=[_primer_pair()],
        max_primer_hits=1,
        max_primer_pair_hits=1,
    )


def test_and_pick_primer_pairs(
    picking_ref: Path,
    zero_score_filtering_params: FilteringParams,
) -> None:
    target = _target()
    pp = _primer_pair(target=target, amplicon_tm=92.96746617835336)
    params = replace(
        zero_score_filtering_params,
        amplicon_tms=MinOptMax(pp.amplicon_tm, pp.amplicon_tm, pp.amplicon_tm),
    )

    offtarget_detector = OffTargetDetector(
        ref=picking_ref,
        max_primer_hits=2,
        max_primer_pair_hits=2,
        three_prime_region_length=5,
        max_mismatches_in_three_prime_region=0,
        max_mismatches=0,
        max_amplicon_size=params.amplicon_sizes.max,
    )

    with pysam.FastaFile(f"{picking_ref}") as fasta:
        designed_primer_pairs = build_and_pick_primer_pairs(
            lefts=[pp.left_primer],
            rights=[pp.right_primer],
            target=target,
            num_primers=1,
            min_difference=1,
            params=params,
            offtarget_detector=offtarget_detector,
            dimer_checker=NtThermoAlign(),
            fasta=fasta,
        )
        assert len(designed_primer_pairs) == 1
        designed_primer_pair = designed_primer_pairs[0]
        assert designed_primer_pair.left_primer == pp.left_primer
        assert designed_primer_pair.right_primer == pp.right_primer
        assert designed_primer_pair.amplicon == pp.amplicon
        assert designed_primer_pair.penalty == pp.penalty
        assert designed_primer_pair.amplicon_tm == pp.amplicon_tm
        assert len(designed_primer_pair.amplicon_sequence) == pp.amplicon.length
