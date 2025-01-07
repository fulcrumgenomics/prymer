import re
from dataclasses import dataclass
from pathlib import Path

import pytest
from fgpyo.sam import Cigar

from prymer.api.oligo import Oligo
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.api.span import Strand
from prymer.offtarget.bwa import BWA_EXECUTABLE_NAME
from prymer.offtarget.bwa import BwaHit
from prymer.offtarget.bwa import BwaResult
from prymer.offtarget.bwa import Query
from prymer.offtarget.offtarget_detector import OffTargetDetector
from prymer.offtarget.offtarget_detector import OffTargetResult


def _build_detector(
    ref_fasta: Path,
    max_primer_hits: int = 1,
    max_primer_pair_hits: int = 1,
    three_prime_region_length: int = 20,
    max_mismatches_in_three_prime_region: int = 0,
    max_mismatches: int = 0,
    min_amplicon_size: int = 10,
    max_amplicon_size: int = 250,
    cache_results: bool = True,
) -> OffTargetDetector:
    """Builds an `OffTargetDetector` with strict defaults"""
    return OffTargetDetector(
        ref=ref_fasta,
        max_primer_hits=max_primer_hits,
        max_primer_pair_hits=max_primer_pair_hits,
        three_prime_region_length=three_prime_region_length,
        max_mismatches_in_three_prime_region=max_mismatches_in_three_prime_region,
        max_mismatches=max_mismatches,
        min_amplicon_size=min_amplicon_size,
        max_amplicon_size=max_amplicon_size,
        cache_results=cache_results,
        keep_spans=True,
        keep_primer_spans=True,
        executable=BWA_EXECUTABLE_NAME,
    )


@pytest.fixture
def multimap_primer_pair() -> PrimerPair:
    """A primer pair that maps to many locations (204 for each primer, 856 as a pair)"""
    return PrimerPair(
        left_primer=Oligo(
            bases="AAAAA",
            tm=37,
            penalty=0,
            span=Span("chr1", start=67, end=71),
        ),
        right_primer=Oligo(
            bases="TTTTT",
            tm=37,
            penalty=0,
            span=Span("chr1", start=75, end=79, strand=Strand.NEGATIVE),
        ),
        amplicon_sequence="AAAAAtacAAAAA",
        amplicon_tm=37,
        penalty=0.0,
    )


# Test that using the cache (or not) does not affect the results
@pytest.mark.parametrize("cache_results", [True, False])
def test_filter(ref_fasta: Path, multimap_primer_pair: PrimerPair, cache_results: bool) -> None:
    primers = list(multimap_primer_pair)
    # keep all, # of mappings found is at the limit
    with _build_detector(
        ref_fasta=ref_fasta, max_primer_hits=204, cache_results=cache_results
    ) as d:
        assert d.filter(primers=primers) == primers
    # keep all, # of mappings found is one beyond the limit
    with _build_detector(
        ref_fasta=ref_fasta, max_primer_hits=203, cache_results=cache_results
    ) as d:
        assert len(d.filter(primers=primers)) == 0


# Test that using the cache (or not) does not affect the results
@pytest.mark.parametrize("cache_results", [True, False])
def test_ok(ref_fasta: Path, multimap_primer_pair: PrimerPair, cache_results: bool) -> None:
    result: OffTargetResult
    with _build_detector(
        ref_fasta=ref_fasta,
        max_primer_hits=204,
        max_primer_pair_hits=856,
        cache_results=cache_results,
    ) as d:
        result = d.check_one(primer_pair=multimap_primer_pair)
        assert result.primer_pair == multimap_primer_pair, id
        assert len(result.left_primer_spans) == 204
        assert len(result.right_primer_spans) == 204
        assert len(result.spans) == 856


# Test that using the cache (or not) does not affect the results
@pytest.mark.parametrize("cache_results", [True, False])
@pytest.mark.parametrize(
    "id, max_primer_hits, max_primer_pair_hits, passes",
    [
        ("too many primer hits", 203, 856, False),
        ("too many primer pair hits", 204, 855, False),
        ("at the maximum", 204, 856, True),
    ],
)
def test_check_too_many_primer_pair_hits(
    ref_fasta: Path,
    multimap_primer_pair: PrimerPair,
    id: str,
    max_primer_hits: int,
    max_primer_pair_hits: int,
    passes: bool,
    cache_results: bool,
) -> None:
    result: OffTargetResult
    with _build_detector(
        ref_fasta=ref_fasta,
        max_primer_hits=max_primer_hits,
        max_primer_pair_hits=max_primer_pair_hits,
        cache_results=cache_results,
    ) as d:
        # when checking caching, calls check_one twice, with the second time the result being.  Do
        # not run twice when we aren't caching, since the second time we map a read, we get
        # different multi-mappings due to how BWA uses a random seed at the start of its execution.
        num_rounds = 2 if cache_results else 1
        for i in range(num_rounds):
            result = d.check_one(primer_pair=multimap_primer_pair)
            assert result.primer_pair == multimap_primer_pair, id
            assert result.passes is passes, id
            # only retrieved from the cache on the first loop iteration, and if using the cache
            cached = cache_results and (i == 1)
            assert result.cached is cached, id


# Test that using the cache (or not) does not affect the results
@pytest.mark.parametrize("cache_results", [True, False])
def test_mappings_of(ref_fasta: Path, cache_results: bool) -> None:
    with _build_detector(ref_fasta=ref_fasta, cache_results=cache_results) as detector:
        p1: Oligo = Oligo(
            tm=37,
            penalty=0,
            span=Span(refname="chr1", start=1, end=30),
            bases="CAGGTGGATCATGAGGTCAGGAGTTCAAGA",
        )
        # NB: the expected hit is returned on the _opposite_ strand
        expected_hit1: BwaHit = BwaHit(
            refname="chr1", start=1, negative=False, cigar=Cigar.from_cigarstring("30M"), edits=0
        )

        p2: Oligo = Oligo(
            tm=37,
            penalty=0,
            span=Span(refname="chr1", start=61, end=93, strand=Strand.NEGATIVE),
            bases="CATGCCCAGCTAATTTTTTGTATTTTTAGTAGA",
        )
        # NB: the expected hit is returned on the _opposite_ strand
        expected_hit2: BwaHit = BwaHit(
            refname="chr1", start=61, negative=True, cigar=Cigar.from_cigarstring("33M"), edits=0
        )

        # Test running the same primers through mappings_of to ensure we get the same results
        for _ in range(10):
            results_dict: dict[str, BwaResult] = detector.mappings_of(primers=[p1, p2])
            assert len(results_dict) == 2
            assert results_dict[p1.bases].hit_count == 1
            assert results_dict[p1.bases].hits[0] == expected_hit1
            assert results_dict[p2.bases].hit_count == 1
            assert results_dict[p2.bases].hits[0] == expected_hit2


# Test building an OffTargetResult for a primer pair with left/right hits on different references
# and in different orientations
def test_build_off_target_result(ref_fasta: Path) -> None:
    hits_by_primer: dict[str, BwaResult] = {
        "A" * 100: BwaResult(
            query=Query(
                id="left",
                bases="A" * 100,
            ),
            hit_count=3,
            hits=[
                BwaHit.build("chr1", 100, False, "100M", 0),
                BwaHit.build("chr1", 400, True, "100M", 0),
                BwaHit.build("chr2", 100, False, "100M", 0),
                BwaHit.build("chr3", 700, True, "100M", 0),
            ],
        ),
        "C" * 100: BwaResult(
            query=Query(
                id="right",
                bases="C" * 100,
            ),
            hit_count=2,
            hits=[
                BwaHit.build("chr1", 800, False, "100M", 0),
                BwaHit.build("chr1", 200, True, "100M", 0),
                BwaHit.build("chr3", 600, False, "100M", 0),
            ],
        ),
    }

    primer_pair = PrimerPair(
        left_primer=Oligo(
            tm=50,
            penalty=0,
            span=Span(refname="chr10", start=100, end=199, strand=Strand.POSITIVE),
            bases="A" * 100,
        ),
        right_primer=Oligo(
            tm=50,
            penalty=0,
            span=Span(refname="chr10", start=300, end=399, strand=Strand.NEGATIVE),
            bases="C" * 100,
        ),
        amplicon_tm=100,
        penalty=0,
    )

    with _build_detector(
        ref_fasta=ref_fasta, max_primer_hits=10, max_primer_pair_hits=10
    ) as detector:
        off_target_result: OffTargetResult = detector._build_off_target_result(
            primer_pair=primer_pair,
            hits_by_primer=hits_by_primer,
        )

    assert set(off_target_result.spans) == {
        Span(refname="chr1", start=100, end=299, strand=Strand.POSITIVE),
        Span(refname="chr3", start=600, end=799, strand=Strand.NEGATIVE),
    }


# Test that using the cache (or not) does not affect the results
@pytest.mark.parametrize("cache_results", [True, False])
@pytest.mark.parametrize(
    "test_id, positive, negative, strand, expected",
    [
        (
            "No mappings - overlapping primers (1bp overlap)",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 199, True, "100M", 0),
            Strand.POSITIVE,
            [],
        ),
        (
            "No mappings - amplicon size too big (1bp too big)",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 151, True, "100M", 0),
            Strand.POSITIVE,
            [],
        ),
        (
            "Mappings - FR pair (R1 F)",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 200, True, "100M", 0),
            Strand.POSITIVE,
            [Span(refname="chr1", start=100, end=299, strand=Strand.POSITIVE)],
        ),
        (
            "Mappings - FR pair (R1 R)",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 200, True, "100M", 0),
            Strand.NEGATIVE,
            [Span(refname="chr1", start=100, end=299, strand=Strand.NEGATIVE)],
        ),
    ],
)
def test_to_amplicons(
    ref_fasta: Path,
    test_id: str,
    positive: BwaHit,
    negative: BwaHit,
    strand: Strand,
    expected: list[Span],
    cache_results: bool,
) -> None:
    with _build_detector(ref_fasta=ref_fasta, cache_results=cache_results) as detector:
        actual = detector._to_amplicons(
            positive_hits=[positive],
            negative_hits=[negative],
            min_len=200,
            max_len=250,
            strand=strand,
        )
        assert actual == expected, test_id


# Test that using the cache (or not) does not affect the results
# NB: BwaHit and Span coordinates are 1-based end-inclusive
#
# One mapping - Overlapping primers (1bp overlap)
# >>>>>>>>>>
#          <<<<<<<<<<
# 19bp amplicon [1,19]
#
# One mapping - Overlapping primers (1bp amplicon)
#          >>>>>>>>>>
# <<<<<<<<<<
# 1bp amplicon [10,10]
#
# No mappings - amplicon length would be 0.
#           >>>>>>>>>>
# <<<<<<<<<<
@pytest.mark.parametrize("cache_results", [True, False])
@pytest.mark.parametrize(
    "test_id, positive, negative, strand, expected",
    [
        (
            "One mapping - Overlapping primers (1bp overlap)",
            BwaHit.build("chr1", 1, False, "10M", 0),
            BwaHit.build("chr1", 10, True, "10M", 0),
            Strand.POSITIVE,
            [Span(refname="chr1", start=1, end=19, strand=Strand.POSITIVE)],
        ),
        (
            "One mapping - Overlapping primers (1bp amplicon)",
            BwaHit.build("chr1", 10, False, "10M", 0),
            BwaHit.build("chr1", 1, True, "10M", 0),
            Strand.POSITIVE,
            [Span(refname="chr1", start=10, end=10, strand=Strand.POSITIVE)],
        ),
        (
            "No mappings",
            BwaHit.build("chr1", 11, False, "10M", 0),
            BwaHit.build("chr1", 1, True, "10M", 0),
            Strand.POSITIVE,
            [],
        ),
    ],
)
def test_to_amplicons_overlapping(
    ref_fasta: Path,
    test_id: str,
    positive: BwaHit,
    negative: BwaHit,
    strand: Strand,
    expected: list[Span],
    cache_results: bool,
) -> None:
    with _build_detector(ref_fasta=ref_fasta, cache_results=cache_results) as detector:
        actual = detector._to_amplicons(
            positive_hits=[positive],
            negative_hits=[negative],
            min_len=1,
            max_len=250,
            strand=strand,
        )
        assert actual == expected, test_id


@pytest.mark.parametrize("cache_results", [True, False])
@pytest.mark.parametrize(
    "positive, negative, expected_error",
    [
        (
            # No mappings - different refnames
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr2", 100, True, "100M", 0),
            "Hits are present on more than one reference",
        ),
        (
            # No mappings - FF pair
            BwaHit.build("chr1", 100, True, "100M", 0),
            BwaHit.build("chr1", 100, True, "100M", 0),
            "Positive hits must be on the positive strand",
        ),
        (
            # No mappings - RR pair
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 100, False, "100M", 0),
            "Negative hits must be on the negative strand",
        ),
    ],
)
def test_to_amplicons_value_error(
    ref_fasta: Path,
    positive: BwaHit,
    negative: BwaHit,
    expected_error: str,
    cache_results: bool,
) -> None:
    with (
        _build_detector(ref_fasta=ref_fasta, cache_results=cache_results) as detector,
        pytest.raises(ValueError, match=expected_error),
    ):
        detector._to_amplicons(
            positive_hits=[positive],
            negative_hits=[negative],
            min_len=200,
            max_len=250,
            strand=Strand.POSITIVE,
        )


def test_generic_filter(ref_fasta: Path) -> None:
    """
    This test isn't intended to validate any runtime assertions, but is a minimal example for the
    type checker to ensure that we can apply `OffTargetDetector.filter()` to arbitrary subclases of
    `Primer`.
    """

    @dataclass(frozen=True)
    class CustomPrimer(Oligo):
        foo: str = "foo"

    # fmt: off
    primers: list[CustomPrimer] = [
        CustomPrimer(bases="AAAA", tm=37, penalty=0, span=Span(refname="chr1", start=1, end=4), foo="bar"),  # noqa: E501
        CustomPrimer(bases="TTTT", tm=37, penalty=0, span=Span(refname="chr2", start=1, end=4), foo="qux"),  # noqa: E501
    ]
    # fmt: on

    with _build_detector(ref_fasta=ref_fasta) as detector:
        # Here we are validating that we can both
        # 1. Pass a list of a `Primer` subclass to `filter()` and
        # 2. Return a list of the same type.
        # NB: we're ignoring the unused value error because we want to check the type hint
        filtered_primers: list[CustomPrimer] = detector.filter(primers)  # noqa: F841


# fmt: off
@pytest.mark.parametrize(
    (
        "max_primer_hits,max_primer_pair_hits,min_primer_pair_hits,three_prime_region_length,"
        "max_mismatches_in_three_prime_region,max_mismatches,max_amplicon_size,min_amplicon_size,"
        "max_gap_opens,max_gap_extends,expected_error"
    ),
    [
        (-1, 1, 1, 20, 0, 0, 1, 1, 0, 0, "'max_primer_hits' must be greater than or equal to 0. Saw -1"),  # noqa: E501
        (1, -1, 1, 20, 0, 0, 1, 1, 0, 0, "'max_primer_pair_hits' must be greater than or equal to 0. Saw -1"),  # noqa: E501
        (1, 1, -1, 20, 0, 0, 1, 1, 0, 0, "'min_primer_pair_hits' must be greater than or equal to 0. Saw -1"),  # noqa: E501
        (1, 1, 1, 5, 0, 0, 1, 1, 0, 0, "'three_prime_region_length' must be greater than or equal to 8. Saw 5"),  # noqa: E501
        (1, 1, 1, 20, -1, 0, 1, 1, 0, 0, "'max_mismatches_in_three_prime_region' must be between 0 and 'three_prime_region_length'=20 inclusive. Saw -1"),  # noqa: E501
        (1, 1, 1, 20, 21, 0, 1, 1, 0, 0, "'max_mismatches_in_three_prime_region' must be between 0 and 'three_prime_region_length'=20 inclusive. Saw 21"),  # noqa: E501
        (1, 1, 1, 20, 0, -1, 1, 1, 0, 0, "'max_mismatches' must be greater than or equal to 0. Saw -1"),  # noqa: E501
        (1, 1, 1, 20, 0, 0, 0, 1, 0, 0, "'max_amplicon_size' must be greater than 0. Saw 0"),
        (1, 1, 1, 20, 0, 0, 1, 1, -1, 0, "'max_gap_opens' must be greater than or equal to 0. Saw -1"), # noqa: E501
        (1, 1, 1, 20, 0, 5, 1, 1, 0, -2, re.escape("'max_gap_extends' must be -1 (for unlimited extensions up to 'max_mismatches'=5) or greater than or equal to 0. Saw -2")), #noqa: E501
        (1, 1, 1, 20, 0, 0, 10, 0, 0, 0, "'min_amplicon_size' must be between 1 and 'max_amplicon_size'=10 inclusive. Saw 0"),  # noqa: E501
        (1, 1, 1, 20, 0, 0, 10, 11, 0, 0, "'min_amplicon_size' must be between 1 and 'max_amplicon_size'=10 inclusive. Saw 11"),  # noqa: E501
    ],
)
# fmt: on
def test_init(
    ref_fasta: Path,
    max_primer_hits: int,
    max_primer_pair_hits: int,
    min_primer_pair_hits: int,
    three_prime_region_length: int,
    max_mismatches_in_three_prime_region: int,
    max_mismatches: int,
    max_amplicon_size: int,
    min_amplicon_size: int,
    max_gap_opens: int,
    max_gap_extends: int,
    expected_error: str,
) -> None:
    with pytest.raises(ValueError, match=expected_error):
        OffTargetDetector(
            ref=ref_fasta,
            max_primer_hits=max_primer_hits,
            max_primer_pair_hits=max_primer_pair_hits,
            min_primer_pair_hits=min_primer_pair_hits,
            three_prime_region_length=three_prime_region_length,
            max_mismatches_in_three_prime_region=max_mismatches_in_three_prime_region,
            max_mismatches=max_mismatches,
            max_amplicon_size=max_amplicon_size,
            min_amplicon_size=min_amplicon_size,
            max_gap_opens=max_gap_opens,
            max_gap_extends=max_gap_extends,
        )
