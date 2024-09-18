from dataclasses import dataclass
from pathlib import Path

import pytest
from fgpyo.sam import Cigar

from prymer.api.primer import Primer
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.api.span import Strand
from prymer.offtarget.bwa import BwaHit
from prymer.offtarget.bwa import BwaResult
from prymer.offtarget.offtarget_detector import OffTargetDetector
from prymer.offtarget.offtarget_detector import OffTargetResult


def _build_detector(
    ref_fasta: Path,
    max_primer_hits: int = 1,
    max_primer_pair_hits: int = 1,
    three_prime_region_length: int = 20,
    max_mismatches_in_three_prime_region: int = 0,
    max_mismatches: int = 0,
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
        max_amplicon_size=max_amplicon_size,
        cache_results=cache_results,
        keep_spans=True,
        keep_primer_spans=True,
    )


@pytest.fixture
def multimap_primer_pair() -> PrimerPair:
    """A primer pair that maps to many locations (204 for each primer, 856 as a pair)"""
    return PrimerPair(
        left_primer=Primer(
            bases="AAAAA",
            tm=37,
            penalty=0,
            span=Span("chr1", start=67, end=71),
        ),
        right_primer=Primer(
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
        p1: Primer = Primer(
            tm=37,
            penalty=0,
            span=Span(refname="chr1", start=1, end=30),
            bases="CAGGTGGATCATGAGGTCAGGAGTTCAAGA",
        )
        # NB: the expected hit is returned on the _opposite_ strand
        expected_hit1: BwaHit = BwaHit(
            refname="chr1", start=1, negative=False, cigar=Cigar.from_cigarstring("30M"), edits=0
        )

        p2: Primer = Primer(
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


# Test that using the cache (or not) does not affect the results
@pytest.mark.parametrize("cache_results", [True, False])
@pytest.mark.parametrize(
    "test_id, left, right, expected",
    [
        (
            "No mappings - different refnames",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr2", 100, True, "100M", 0),
            [],
        ),
        (
            "No mappings - FF pair",
            BwaHit.build("chr1", 100, True, "100M", 0),
            BwaHit.build("chr1", 100, True, "100M", 0),
            [],
        ),
        (
            "No mappings - RR pair",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 100, False, "100M", 0),
            [],
        ),
        (
            "No mappings - overlapping primers (1bp overlap)",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 199, True, "100M", 0),
            [],
        ),
        (
            "No mappings - amplicon size too big (1bp too big)",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 151, True, "100M", 0),
            [],
        ),
        (
            "Mappings - FR pair (R1 F)",
            BwaHit.build("chr1", 100, False, "100M", 0),
            BwaHit.build("chr1", 200, True, "100M", 0),
            [Span(refname="chr1", start=100, end=299)],
        ),
        (
            "Mappings - FR pair (R1 R)",
            BwaHit.build("chr1", 200, True, "100M", 0),
            BwaHit.build("chr1", 100, False, "100M", 0),
            [Span(refname="chr1", start=100, end=299)],
        ),
    ],
)
def test_to_amplicons(
    ref_fasta: Path,
    test_id: str,
    left: BwaHit,
    right: BwaHit,
    expected: list[Span],
    cache_results: bool,
) -> None:
    with _build_detector(ref_fasta=ref_fasta, cache_results=cache_results) as detector:
        actual = detector._to_amplicons(lefts=[left], rights=[right], max_len=250)
        assert actual == expected, test_id


def test_generic_filter(ref_fasta: Path) -> None:
    """
    This test isn't intended to validate any runtime assertions, but is a minimal example for the
    type checker to ensure that we can apply `OffTargetDetector.filter()` to arbitrary subclases of
    `Primer`.
    """

    @dataclass(frozen=True)
    class CustomPrimer(Primer):
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
