import pytest
from fgpyo.fasta.sequence_dictionary import SequenceDictionary

from prymer.model import Span
from prymer.model import Strand


@pytest.mark.parametrize(
    "refname, invalid_start, invalid_end, strand",
    [
        ("chr1", 0, 10, "+"),  # start < 1
        ("chr1", 30, 20, "-"),  # start > end
        (" ", 0, 10, "+"),  # empty refname
    ],
)
def test_invalid_span(refname: str, invalid_start: int, invalid_end: int, strand: str) -> None:
    """Test that invalid spans raise an error"""
    with pytest.raises(ValueError):
        Span(refname, invalid_start, invalid_end, Strand(strand))


@pytest.mark.parametrize(
    "base_idx, expected_rel_coords",
    [
        (100, 0),  # 1 based start
        (101, 1),  # second base
        (150, 50),  # span midpoint
        (199, 99),  # second to last base
        (200, 100),  # 1 based + inclusive end
    ],
)
def test_get_offset_valid(base_idx: int, expected_rel_coords: int) -> None:
    """Test whether a Span position is correctly converted to a relative coord"""
    span = Span("chr1", 100, 200, Strand.POSITIVE)
    assert span.get_offset(position=base_idx) == expected_rel_coords


@pytest.mark.parametrize(
    "base_idx",
    [
        (999),  # base_idx < self.start
        (2001),  # base_idx > self.end
    ],
)
def test_get_offset_raises(base_idx: int) -> None:
    """Test whether get_offset() raises error with invalid inputs"""
    span = Span("chr1", 1000, 2000, Strand.POSITIVE)
    with pytest.raises(ValueError):
        span.get_offset(position=base_idx)


@pytest.mark.parametrize(
    "offset, submap_length, expected_start, expected_end",
    [
        (0, 1, 1, 1),  # 1-based inclusive, start == end
        (1, 2, 2, 3),  # wholly contained
        (5, 5, 6, 10),  # 1-based start + offset of 6 = 10 == length
    ],
)
def test_get_submap_from_offset_valid(
    offset: int, submap_length: int, expected_start: int, expected_end: int
) -> None:
    """Test whether relative sub-span coords are correctly converted to genomic coords"""

    test_map = Span(refname="chr1", start=1, end=10, strand=Strand.POSITIVE)
    assert test_map.length == 10
    assert test_map.get_subspan(offset=offset, subspan_length=submap_length) == Span(
        refname="chr1", start=expected_start, end=expected_end, strand=Strand.POSITIVE
    )


@pytest.mark.parametrize(
    "offset, submap_length",
    [
        (-1, 2),  # offset < 0
        (5, 0),  # length < 1
        (12, 1),  # offset > length
        (2, 11),  # length outside source span
        (6, 6),  # end of sub-span outside source span
    ],
)
def test_get_subspan_from_offset_raises(offset: int, submap_length: int) -> None:
    """Test if get_submap_from_offset() raises error with invalid sub-span coords"""

    test_map = Span(refname="chr1", start=1, end=10, strand=Strand.POSITIVE)
    with pytest.raises(ValueError):
        test_map.get_subspan(offset=offset, subspan_length=submap_length)


def test_round_trip_conversion() -> None:
    """Test if span <-> subspan modifies positional data unexpectedly"""
    span = Span(refname="chr1", start=100, end=200, strand=Strand.POSITIVE)
    offset = span.get_offset(position=span.start)
    length = span.length
    test_submap = span.get_subspan(offset=offset, subspan_length=length)
    assert span == test_submap


@pytest.mark.parametrize(
    "line, expected_span, expected_length",
    [
        ("chr1:1-10:+", Span(refname="chr1", start=1, end=10), 10),
        ("chr1:1-10", Span(refname="chr1", start=1, end=10), 10),
        ("chr2:20-30:-", Span(refname="chr2", start=20, end=30, strand=Strand.NEGATIVE), 11),
        ("chr2:20 - 30", Span(refname="chr2", start=20, end=30), 11),
        ("chr1:1-10:--", Span(refname="chr1", start=1, end=10, strand=Strand.NEGATIVE), 10),
    ],
)
def test_span_from_valid_string(line: str, expected_span: Span, expected_length: int) -> None:
    """Test whether from_string() yields expected Span objects"""
    span = Span.from_string(line)
    assert span == expected_span
    assert span.length == expected_length


@pytest.mark.parametrize(
    "invalid_line",
    [
        ("chr1:1-10:+:abc"),
        ("chr2:20-30:x"),  # invalid strand
        ("chr1:1-10:--:abc:def"),  # >3 colon-delimited fields
        ("chr1:1-x:+"),  # string end
    ],
)
def test_spanfrom_invalid_string_raises(invalid_line: str) -> None:
    """Test whether from_string() yields Span objects as expected from string input"""
    with pytest.raises((ValueError, TypeError)):
        Span.from_string(invalid_line)


@pytest.mark.parametrize(
    "span1, span2, length_of_overlap",
    [
        (Span("chr1", 500, 1000), Span("chr1", 500, 1000), 501),
        (Span("chr1", 500, 1000), Span("chr1", 600, 700), 101),
        (Span("chr1", 500, 1000), Span("chr1", 400, 2000), 501),
        (Span("chr1", 500, 1000), Span("chr1", 750, 850), 101),
        (Span("chr1", 500, 1000), Span("chr1", 1000, 1001), 1),
        (Span("chr1", 500, 1000), Span("chr2", 500, 1000), 0),
        (Span("chr1", 500, 1000), Span("chr1", 5000, 10000), 0),
        (Span("chr1", 500, 1000), Span("chr1", 1, 499), 0),
        (
            Span("chr1", 500, 1000),
            Span("chr1", 1001, 1100, strand=Strand.NEGATIVE),
            0,
        ),
    ],
)
def test_span_overlap(span1: Span, span2: Span, length_of_overlap: int) -> None:
    overlaps: bool = length_of_overlap > 0
    assert span1.overlaps(span2) == overlaps
    assert span2.overlaps(span1) == overlaps
    assert span1.length_of_overlap_with(span2) == length_of_overlap
    assert span2.length_of_overlap_with(span1) == length_of_overlap


@pytest.mark.parametrize(
    "this, that, expected",
    [
        # same span
        (Span("chr1", 100, 200), Span("chr1", 100, 200), 0),
        (Span("chr1", 100, 200, Strand.POSITIVE), Span("chr1", 100, 200, Strand.POSITIVE), 0),
        (Span("chr1", 100, 200, Strand.NEGATIVE), Span("chr1", 100, 200, Strand.NEGATIVE), 0),
        # earlier reference
        (Span("chr1", 100, 200, Strand.NEGATIVE), Span("chr2", 1, 2, Strand.POSITIVE), -1),
        (Span("chr1", 100, 200, Strand.NEGATIVE), Span("chr3", 1, 2, Strand.POSITIVE), -1),
        (Span("chr2", 100, 200, Strand.NEGATIVE), Span("chr3", 1, 2, Strand.POSITIVE), -1),
        # same reference, earlier start
        (Span("chr1", 99, 200, Strand.NEGATIVE), Span("chr1", 100, 150, Strand.POSITIVE), -1),
        # same reference and start, earlier end
        (Span("chr1", 100, 199, Strand.NEGATIVE), Span("chr1", 100, 200, Strand.POSITIVE), -1),
        # same reference, start, and end, but earlier strand
        (Span("chr1", 100, 200, Strand.POSITIVE), Span("chr1", 100, 200, Strand.NEGATIVE), -1),
    ],
)
def test_span_compare(this: Span, that: Span, expected: int, seq_dict: SequenceDictionary) -> None:
    assert expected == Span.compare(this=this, that=that, seq_dict=seq_dict)
    assert -expected == Span.compare(this=that, that=this, seq_dict=seq_dict)
