from dataclasses import dataclass
from dataclasses import replace
from typing import Optional

import pytest
from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.sequence import reverse_complement

from prymer.api.oligo import Oligo
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.api.span import Strand


@dataclass(init=True, frozen=True)
class PrimerPairTestCase:
    """Test case for a `PrimerPair`.

    Attributes:
        primer_pair: the PrimerPair to test
        bed_12_fields: the fields, that when tab-delimited, are the expected string
            for the `PrimerPair.bed_12` method
        gc_pct: the expected value for the `PrimerPair.percent_gc_content` property
        inner: the expected value for the `PrimerPair.inner` method
        length: the expected value for the `PrimerPair.length` property
        str_fields: the fields, that when tab-delimited, are the expected string for the
            `Primer.__str__` method
    """

    primer_pair: PrimerPair
    gc_pct: float
    length: int

    @staticmethod
    def primer_pair_from_left_primer(left: Oligo, right_offset: int = 50) -> PrimerPair:
        """
        Generates a PrimerPair for use in unit tests. Will first generate
        a right primer using a standard formula, and then calculates the
        PrimerPair fields based on the left & right primers.
        """
        right: Oligo = PrimerPairTestCase.right_primer_from_left_primer(
            left=left, right_offset=right_offset
        )

        amplicon = replace(left.span, end=right.span.end)
        amplicon_sequence = PrimerPairTestCase._to_amplicon_sequence(
            left.bases,
            right.bases,
            amplicon.length - 1,
        )
        amplicon_tm = left.tm + right.tm
        penalty = left.penalty + right.penalty

        return PrimerPair(
            left_primer=left,
            right_primer=right,
            amplicon_sequence=amplicon_sequence,
            amplicon_tm=amplicon_tm,
            penalty=penalty,
        )

    @staticmethod
    def right_primer_from_left_primer(left: Oligo, right_offset: int) -> Oligo:
        """
        Provides a standard conversion for a left primer to a right primer for use in
        tests of PrimerPair.
        """
        return replace(
            left,
            bases=reverse_complement(left.bases) if left.bases is not None else None,
            tm=left.tm + 10,
            penalty=left.penalty + 10,
            span=replace(
                left.span,
                start=left.span.start + right_offset,
                end=left.span.end + right_offset,
                strand=Strand.NEGATIVE,
            ),
        )

    @staticmethod
    def _to_amplicon_sequence(
        left: Optional[str], right: Optional[str], length: int
    ) -> Optional[str]:
        """
        Provides a consistent mechanism to generate a valid amplicon sequence.
        There's no inherent meaning to the algorithm.
        """
        if left is None or right is None:
            return None

        # Combine the bases from left and right primers
        bases = "".join(
            f"{left_base}{right_base}" for left_base, right_base in zip(left, right, strict=True)
        )
        # Repeat the combined bases enough times to cover the required length
        times = (length // len(bases)) + 1
        extended_bases = (bases * times)[: length + 1]
        return extended_bases


def build_primer_pair_test_cases() -> list[PrimerPairTestCase]:
    """
    Builds a set of test cases for `PrimerPair` methods using a
    supplied left primer.
    """
    return [
        PrimerPairTestCase(
            primer_pair=PrimerPairTestCase.primer_pair_from_left_primer(
                left=Oligo(
                    bases="GATTACA",
                    tm=12.34,
                    penalty=56.78,
                    span=Span(refname="chr1", start=1, end=7, strand=Strand.POSITIVE),
                )
            ),
            gc_pct=29.825,
            length=57,
        ),
        PrimerPairTestCase(
            primer_pair=PrimerPairTestCase.primer_pair_from_left_primer(
                left=Oligo(
                    bases="TGTAATC",
                    tm=87.65,
                    penalty=43.21,
                    span=Span(refname="chr22", start=100, end=106, strand=Strand.POSITIVE),
                )
            ),
            gc_pct=28.07,
            length=57,
        ),
        PrimerPairTestCase(
            primer_pair=PrimerPairTestCase.primer_pair_from_left_primer(
                left=Oligo(
                    bases="A" * 40,
                    tm=12.34,
                    penalty=56.78,
                    span=Span(refname="chr1", start=1, end=40, strand=Strand.POSITIVE),
                )
            ),
            gc_pct=0.0,
            length=90,
        ),
        PrimerPairTestCase(
            primer_pair=PrimerPairTestCase.primer_pair_from_left_primer(
                left=Oligo(
                    bases="GGGGGGG",
                    tm=12.34,
                    penalty=56.78,
                    span=Span(refname="chr1", start=1, end=7, strand=Strand.POSITIVE),
                )
            ),
            gc_pct=100.0,
            length=57,
        ),
        # This one is written out long form in order to override the automated
        # generation of the right primer. We want to handle a case where the primers
        # overlap
        PrimerPairTestCase(
            primer_pair=PrimerPair(
                left_primer=Oligo(
                    bases="GATTACA",
                    tm=12.34,
                    penalty=56.78,
                    span=Span(refname="chr1", start=1, end=7, strand=Strand.POSITIVE),
                ),
                right_primer=Oligo(
                    bases="TGTAATC",
                    tm=87.65,
                    penalty=43.21,
                    span=Span(refname="chr1", start=5, end=11, strand=Strand.NEGATIVE),
                ),
                amplicon_tm=25,
                penalty=99.99,
            ),
            gc_pct=0,
            length=11,
        ),
    ]


PRIMER_PAIR_TEST_CASES = build_primer_pair_test_cases()


@pytest.mark.parametrize("test_case", PRIMER_PAIR_TEST_CASES)
def test_primer_pair_gc_content(test_case: PrimerPairTestCase) -> None:
    """Test the `PrimerPair.percent_gc_content` property."""
    assert test_case.primer_pair.percent_gc_content == test_case.gc_pct


@pytest.mark.parametrize("test_case", PRIMER_PAIR_TEST_CASES)
def test_primer_pair_length(test_case: PrimerPairTestCase) -> None:
    """Test the `PrimerPair.length` property"""
    assert test_case.primer_pair.length == test_case.length


@pytest.mark.parametrize("test_case", PRIMER_PAIR_TEST_CASES)
def test_primer_pair_span(test_case: PrimerPairTestCase) -> None:
    """
    Test that the `PrimerPair.span property returns the amplicon field
    """
    assert test_case.primer_pair.span == test_case.primer_pair.amplicon


@pytest.mark.parametrize("test_case", PRIMER_PAIR_TEST_CASES)
def test_primer_pair_bases(test_case: PrimerPairTestCase) -> None:
    """
    Test that the `PrimerPair.bases property returns the amplicon_sequence field
    """
    assert test_case.primer_pair.bases == test_case.primer_pair.amplicon_sequence


@pytest.mark.parametrize(
    "left_tail,right_tail",
    [
        ("AAA", "GGG"),
        ("GGG", "AAAAA"),
    ],
)
def test_with_tails(left_tail: str, right_tail: str) -> None:
    """
    Test that the with_tails method adds the correct tail(s) to
    the PrimerPair and both Primers
    """
    pp = PRIMER_PAIR_TEST_CASES[0].primer_pair

    assert pp.left_primer.tail is None
    assert pp.right_primer.tail is None

    pp_with_tails = pp.with_tails(left_tail, right_tail)

    assert pp_with_tails.left_primer.tail == left_tail
    assert pp_with_tails.right_primer.tail == right_tail

    pp_empty_str = pp_with_tails.with_tails("", "")
    assert pp_empty_str.left_primer.tail == ""
    assert pp_empty_str.right_primer.tail == ""


def test_with_names() -> None:
    """Test that the with_names method adds names to the PrimerPair and both Primers"""
    # create a primer pair with no names
    pp = PRIMER_PAIR_TEST_CASES[0].primer_pair
    assert pp.name is None
    assert pp.left_primer.name is None
    assert pp.right_primer.name is None

    # add names to all of them
    pp_with_names = pp.with_names("pp_name", "lp_name", "rp_name")
    assert pp_with_names.name == "pp_name"
    assert pp_with_names.left_primer.name == "lp_name"
    assert pp_with_names.right_primer.name == "rp_name"

    # set them all to _empty_ names
    pp_reset = pp_with_names.with_names("", "", "")
    assert pp_reset.name == ""
    assert pp_reset.left_primer.name == ""
    assert pp_reset.right_primer.name == ""


def test_reference_mismatch() -> None:
    """
    Test that a PrimerPair and both Primers all have a Span with
    the same reference sequence
    """

    pp = PRIMER_PAIR_TEST_CASES[0].primer_pair

    with pytest.raises(ValueError, match="different references"):
        replace(
            pp,
            left_primer=replace(
                pp.left_primer,
                span=replace(pp.left_primer.span, refname="no-name"),
            ),
        )

    with pytest.raises(ValueError, match="different references"):
        replace(
            pp,
            right_primer=replace(
                pp.right_primer,
                span=replace(pp.right_primer.span, refname="no-name"),
            ),
        )


def test_right_primer_before_left_primer() -> None:
    """Test that an exception is raised if the left primer starts after the right primer ends"""
    pp = PRIMER_PAIR_TEST_CASES[0].primer_pair
    with pytest.raises(ValueError, match="Left primer does not start before the right primer"):
        replace(
            pp,
            left_primer=pp.right_primer,
            right_primer=pp.left_primer,
        )


def test_iter() -> None:
    """Test that the iter dunder method returns the left and right primers"""
    pp = PRIMER_PAIR_TEST_CASES[0].primer_pair
    assert list(iter(pp)) == [pp.left_primer, pp.right_primer]


@pytest.mark.parametrize(
    "this, that, expected_by_amplicon_true, expected_by_amplicon_false",
    [
        # same primer
        (
            PrimerPairTestCase.primer_pair_from_left_primer(
                Oligo(
                    bases="GATTACA",
                    tm=1.0,
                    penalty=2.0,
                    span=Span(refname="chr1", start=100, end=106, strand=Strand.POSITIVE),
                )
            ),
            PrimerPairTestCase.primer_pair_from_left_primer(
                Oligo(
                    bases="GATTACA",
                    tm=1.0,
                    penalty=2.0,
                    span=Span(refname="chr1", start=100, end=106, strand=Strand.POSITIVE),
                )
            ),
            0,
            0,
        ),
        # different primer (chromosome)
        (
            PrimerPairTestCase.primer_pair_from_left_primer(
                Oligo(
                    bases="GATTACA",
                    tm=1.0,
                    penalty=2.0,
                    span=Span(refname="chr1", start=100, end=106, strand=Strand.POSITIVE),
                )
            ),
            PrimerPairTestCase.primer_pair_from_left_primer(
                Oligo(
                    bases="GATTACA",
                    tm=1.0,
                    penalty=2.0,
                    span=Span(refname="chr2", start=100, end=106, strand=Strand.POSITIVE),
                )
            ),
            -1,
            -1,
        ),
        # same primer when by amplicon, but different by primer
        (
            PrimerPairTestCase.primer_pair_from_left_primer(
                Oligo(
                    bases="GATTAC",
                    tm=1.0,
                    penalty=2.0,
                    span=Span(refname="chr1", start=100, end=105, strand=Strand.POSITIVE),
                ),
                right_offset=51,
            ),
            PrimerPairTestCase.primer_pair_from_left_primer(
                Oligo(
                    bases="GATTACA",
                    tm=1.0,
                    penalty=2.0,
                    span=Span(refname="chr1", start=100, end=106, strand=Strand.POSITIVE),
                ),
                right_offset=50,
            ),
            0,
            -1,
        ),
    ],
)
def test_primer_pair_compare(
    this: PrimerPair,
    that: PrimerPair,
    expected_by_amplicon_true: int,
    expected_by_amplicon_false: int,
    seq_dict: SequenceDictionary,
) -> None:
    # expected_by_amplicon_true
    assert expected_by_amplicon_true == PrimerPair.compare(
        this=this, that=that, seq_dict=seq_dict, by_amplicon=True
    )
    assert -expected_by_amplicon_true == PrimerPair.compare(
        this=that, that=this, seq_dict=seq_dict, by_amplicon=True
    )
    # expected_by_amplicon_false
    assert expected_by_amplicon_false == PrimerPair.compare(
        this=this, that=that, seq_dict=seq_dict, by_amplicon=False
    )
    assert -expected_by_amplicon_false == PrimerPair.compare(
        this=that, that=this, seq_dict=seq_dict, by_amplicon=False
    )


def test_calculate_amplicon_span() -> None:
    left = Oligo(name="l", bases="AACCGGTTAA", tm=60, penalty=1, span=Span("chr1", 50, 59))
    right = Oligo(name="l", bases="AACCGGTTAA", tm=60, penalty=1, span=Span("chr1", 150, 159))
    assert PrimerPair.calculate_amplicon_span(left, right) == Span("chr1", 50, 159)

    left = Oligo(name="l", bases="AACCGGTTAA", tm=60, penalty=1, span=Span("chr2", 50, 59))
    right = Oligo(name="l", bases="AACCGGTTAA", tm=60, penalty=1, span=Span("chr3", 150, 159))
    with pytest.raises(ValueError, match="different references"):
        PrimerPair.calculate_amplicon_span(left, right)

    left = Oligo(name="l", bases="AACCGGTTAA", tm=60, penalty=1, span=Span("chr1", 150, 159))
    right = Oligo(name="l", bases="AACCGGTTAA", tm=60, penalty=1, span=Span("chr1", 50, 59))
    with pytest.raises(ValueError, match="Left primer does not start before the right primer"):
        PrimerPair.calculate_amplicon_span(left, right)

    left = Oligo(name="l", bases="AACCGGTTAAACGTT", tm=60, penalty=1, span=Span("chr1", 150, 164))
    right = Oligo(name="l", bases="AACCGGTTAA", tm=60, penalty=1, span=Span("chr1", 150, 159))
    with pytest.raises(ValueError, match="Right primer ends before left primer ends"):
        PrimerPair.calculate_amplicon_span(left, right)
