from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Optional

import pytest
from fgpyo.fasta.sequence_dictionary import SequenceDictionary

from prymer.api.oligo import Oligo
from prymer.api.span import Span
from prymer.api.span import Strand


@pytest.mark.parametrize(
    "bases,tm,penalty,test_span",
    [
        ("AGCT", 1.0, 2.0, Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE)),
        (
            "AGCTAGCTAA",
            10.0,
            20.0,
            Span(refname="chr2", start=1, end=10, strand=Strand.NEGATIVE),
        ),
    ],
)
def test_valid_primer_construction(bases: str, tm: float, penalty: float, test_span: Span) -> None:
    """Test Oligo construction with valid input and ensure reported lengths match"""
    test_primer = Oligo(bases=bases, tm=tm, penalty=penalty, span=test_span)
    assert test_primer.length == test_primer.span.length


def test_span_returns_span(test_span: Span) -> None:
    """Test that the mapping property returns the span object."""
    test_primer = Oligo(
        bases="AGCTAGCTAA",
        tm=1.0,
        penalty=2.0,
        span=test_span,
    )
    assert test_primer.span == test_span


def test_invalid_primer_construction_raises() -> None:
    """Test Oligo construction with invalid input raises ValueError"""
    with pytest.raises(ValueError, match="Bases must not be an empty string"):
        Oligo(
            bases="",
            tm=1.0,
            penalty=2.0,
            span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
        )

    with pytest.raises(
        ValueError, match="Conflicting lengths"
    ):
        Oligo(
            bases="ACGT",
            tm=1.0,
            penalty=2.0,
            span=Span(refname="chr1", start=1, end=1000, strand=Strand.POSITIVE),
        )


@dataclass(init=True, frozen=True)
class OligoTestCase:
    """Test case for an `Oligo`.

    Attributes:
        primer: the primer to test
        gc_pct: the expected value for the `Oligo.percent_gc_content` method
        longest_hp: the expected value for the `Oligo.longest_homopolymer` method
        longest_dinuc: the expected value for the `Oligo.longest_dinucleotide_run` method
    """

    primer: Oligo
    gc_pct: float
    longest_hp: int
    longest_dinuc: int


def build_primer_test_cases() -> list[OligoTestCase]:
    """Builds a set of test cases for `Primer` methods."""
    return [
        OligoTestCase(
            primer=Oligo(
                bases="ATAT",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
            ),
            gc_pct=0.0,
            longest_hp=1,
            longest_dinuc=4,
        ),
        OligoTestCase(
            primer=Oligo(
                bases="ACGTAAAAAATT",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=12, strand=Strand.POSITIVE),
            ),
            gc_pct=16.667,
            longest_hp=6,
            longest_dinuc=6,
        ),
        OligoTestCase(
            primer=Oligo(
                bases="ATAC",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
            ),
            gc_pct=25.0,
            longest_hp=1,
            longest_dinuc=2,
        ),
        OligoTestCase(
            primer=Oligo(
                bases="ATATCC",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=6, strand=Strand.POSITIVE),
            ),
            gc_pct=33.333,
            longest_hp=2,
            longest_dinuc=4,
        ),
        OligoTestCase(
            primer=Oligo(
                bases="AGCT",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
            ),
            gc_pct=50.0,
            longest_hp=1,
            longest_dinuc=2,
        ),
        OligoTestCase(
            primer=Oligo(
                bases="GGGGG",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=5, strand=Strand.POSITIVE),
            ),
            gc_pct=100.0,
            longest_hp=5,
            longest_dinuc=4,
        ),
        OligoTestCase(
            primer=Oligo(
                bases="ccgTATGC",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=8, strand=Strand.POSITIVE),
            ),
            gc_pct=62.5,
            longest_hp=2,
            longest_dinuc=2,
        ),
        OligoTestCase(
            primer=Oligo(
                bases="ACGT",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
            ),
            gc_pct=50.0,
            longest_hp=1,
            longest_dinuc=0,
        ),
        OligoTestCase(
            primer=Oligo(
                bases="ACACACTCTCTCT",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=1, end=13, strand=Strand.POSITIVE),
            ),
            gc_pct=46.154,
            longest_hp=1,
            longest_dinuc=8,
        ),
    ]


OLIGO_TEST_CASES: list[OligoTestCase] = build_primer_test_cases()


@pytest.mark.parametrize("test_case", OLIGO_TEST_CASES)
def test_gc_content_calc(test_case: OligoTestCase) -> None:
    """Test that percent GC content is calculated correctly."""
    assert test_case.primer.percent_gc_content == pytest.approx(test_case.gc_pct)


@pytest.mark.parametrize("test_case", OLIGO_TEST_CASES)
def test_longest_homopolymer_len_calc(test_case: OligoTestCase) -> None:
    """Test that longest homopolymer run is calculated correctly."""
    assert test_case.primer.longest_hp_length == test_case.longest_hp


@pytest.mark.parametrize(
    "init, value, expected",
    [
        # no initial value
        (None, "", ""),
        (None, "GATTACA", "GATTACA"),
        (None, "NNNNN", "NNNNN"),
        # update the initial value
        ("TTTT", "", ""),
        ("TTTT", "GATTACA", "GATTACA"),
        ("TTTT", "NNNNN", "NNNNN"),
    ],
)
def test_with_tail(init: Optional[str], value: str, expected: Optional[str]) -> None:
    """Tests the `with_tail` method, setting the initial value to `init`, updating the
    tail using the `with_tail()` method with value `value`, and testing for the execpted
    value `expected`."""
    test_primer = Oligo(
        bases="AGCT",
        tm=1.0,
        penalty=2.0,
        span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
        tail=init,
    )
    modified_test_primer = test_primer.with_tail(tail=value)
    assert modified_test_primer.tail is expected
    assert modified_test_primer.bases == "AGCT"


@pytest.mark.parametrize(
    "tail_seq, bases, expected_result",
    [
        ("TTTT", "AGCT", "TTTTAGCT"),
        ("", "AGCT", "AGCT"),
        (None, "AGCT", "AGCT"),
        ("NNNNNNNNNN", "AGCT", "NNNNNNNNNNAGCT"),
        ("GATTACA", "AGCT", "GATTACAAGCT"),
    ],
)
def test_bases_with_tail(
    tail_seq: Optional[str], bases: str, expected_result: Optional[str]
) -> None:
    test_primer = Oligo(
        bases=bases,
        tm=1.0,
        penalty=2.0,
        span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
        tail=tail_seq,
    )
    assert test_primer.bases_with_tail == expected_result


@pytest.mark.parametrize(
    "init, value, expected",
    [
        # no initial value
        (None, "", ""),
        (None, "GATTACA", "GATTACA"),
        (None, "NNNNN", "NNNNN"),
        # update the initial value
        ("TTTT", "", ""),
        ("TTTT", "GATTACA", "GATTACA"),
        ("TTTT", "NNNNN", "NNNNN"),
    ],
)
def test_with_name(init: Optional[str], value: str, expected: Optional[str]) -> None:
    test_primer = Oligo(
        bases="AGCT",
        tm=1.0,
        penalty=2.0,
        span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
        name=init,
    )
    modified_test_primer = test_primer.with_name(name=value)
    assert test_primer.name is init
    assert modified_test_primer.name is expected


@pytest.mark.parametrize(
    "name, expected_id",
    [
        ("test", "test"),
        (None, "chr1:1-10:+"),
    ],
)
def test_id_generation(
    test_span: Span,
    name: Optional[str],
    expected_id: str,
) -> None:
    """Asserts that the id field is correctly generated based on the name field."""

    # For each scenario, generate a Primer object and assert that the generated ID
    # matches the expected ID.
    primer = Oligo(
        name=name,
        span=test_span,
        bases="AAAAAAAAAA",
        tm=1.0,
        penalty=1.0,
        tail="AAA",
    )
    assert primer.id == expected_id


@pytest.mark.parametrize("test_case", OLIGO_TEST_CASES)
def test_untailed_length(test_case: OligoTestCase) -> None:
    assert test_case.primer.length == test_case.primer.untailed_length


@pytest.mark.parametrize(
    "tail_seq, expected_length",
    [
        ("TTTT", 8),
        ("", 4),
        (None, 4),
        ("NNNNNNNNNN", 14),
        ("GATTACA", 11),
    ],
)
def test_tailed_length(tail_seq: str, expected_length: int) -> None:
    test_primer = Oligo(
        bases="AGCT",
        tm=1.0,
        penalty=2.0,
        span=Span(refname="chr1", start=1, end=4, strand=Strand.POSITIVE),
        tail=tail_seq,
    )
    assert test_primer.tailed_length == expected_length


def test_primer_serialization_roundtrip() -> None:
    input_primers: list[Oligo] = [test_case.primer for test_case in OLIGO_TEST_CASES]

    with NamedTemporaryFile(suffix=".txt", mode="r", delete=True) as write_file:
        path = Path(write_file.name)

        # write them to a file
        Oligo.write(path, *input_primers)

        # read them back in again
        output_primers = list(Oligo.read(path=path))

        # make sure they're the same!
        assert input_primers == output_primers


@pytest.mark.parametrize(
    "this, that, expected",
    [
        # same primer
        (
            Oligo(
                bases="GATTACA",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=100, end=106, strand=Strand.POSITIVE),
            ),
            Oligo(
                bases="GATTACA",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=100, end=106, strand=Strand.POSITIVE),
            ),
            0,
        ),
        # different primer (chromosome)
        (
            Oligo(
                bases="GATTACA",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr1", start=100, end=106, strand=Strand.POSITIVE),
            ),
            Oligo(
                bases="GATTACA",
                tm=1.0,
                penalty=2.0,
                span=Span(refname="chr2", start=100, end=106, strand=Strand.POSITIVE),
            ),
            -1,
        ),
    ],
)
def test_primer_compare(
    this: Oligo, that: Oligo, expected: int, seq_dict: SequenceDictionary
) -> None:
    assert expected == Oligo.compare(this=this, that=that, seq_dict=seq_dict)
    assert -expected == Oligo.compare(this=that, that=this, seq_dict=seq_dict)
