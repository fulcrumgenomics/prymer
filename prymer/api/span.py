"""
# Span Classes and Methods

The [`Span`][prymer.api.span.Span] class and associated methods represent positional
information from a span of some sequence (like a primer or an amplicon) over a reference sequence
(e.g. genome, target, amplicon).

Class attributes include `refname`, `start`, `end`, and `strand`.

The [`Span`][prymer.api.span.Span] class uses 1-based, closed intervals for start and end
positions.  The [`get_bedlike_coords()`][prymer.api.span.Span.get_bedlike_coords] method
may be used to conver a `Span` to zero-based open-ended coordinates.

## Examples of interacting with the `Span` class


```python
>>> span_string = "chr1:1-10:+"
>>> span = Span.from_string(span_string)
>>> span
Span(refname='chr1', start=1, end=10, strand=<Strand.POSITIVE: '+'>)
>>> print(span.start)
1
>>> # get 0-based position within span at position 10
>>> span.get_offset(position=10)
9
>>> # create a new subspan derived from the source map
>>> new_subspan = span.get_subspan(offset=5, subspan_length=5)
>>> new_subspan
Span(refname='chr1', start=6, end=10, strand=<Strand.POSITIVE: '+'>)
>>> print(span.length)
10

```

A `Span` can be compared to another `Span` to test for overlap, return the length of overlap, and
test if one span fully contains the other:

```python
>>> span1 = Span(refname="chr1", start=50, end=100)
>>> span2 = Span(refname="chr1", start=75, end=90)
>>> span1.overlaps(span2)
True
>>> span2.overlaps(span1)
True
>>> span1.length_of_overlap_with(span2)
16
>>> span1.length_of_overlap_with(Span(refname="chr1", start=75, end=125))
26
>>> span1.overlaps(Span(refname="chr1", start=200, end=225))
False
>>> span1.contains(Span(refname="chr1", start=75, end=125))
False
>>> span1.contains(span2)
True

```

In some cases, it's useful to have the coordinates be converted to zero-based open-ended, for
example with use with the `pysam` module when using the `fetch()` methods for `pysam.AlignmentFile`,
`pysam.FastaFile`, and `pysam.VariantFile`.

```python
>>> span.get_bedlike_coords()
BedLikeCoords(start=0, end=10)

```
"""

from dataclasses import dataclass
from dataclasses import replace
from enum import StrEnum
from enum import unique
from typing import Optional
from typing import Self
from typing import Tuple

from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.util.metric import Metric

from prymer.api import coordmath


@unique
class Strand(StrEnum):
    """Represents the strand of a span to the genome."""

    POSITIVE = "+"
    NEGATIVE = "-"


@dataclass(init=True, frozen=True)
class BedLikeCoords:
    """Represents the coordinates (only, no refname or strand) that correspond to a
    zero-based open-ended position."""

    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"Start position must be >=0; received start={self.start}")
        if self.end < self.start:
            raise ValueError(f"End must be >= start; received start={self.start}, end={self.end}")


@dataclass(init=True, frozen=True, eq=True)
class Span(Metric["Span"]):
    """Represents the span of some sequence (target, primer, amplicon, etc.) to a genome.

    Attributes:
        refname: name of a reference sequence or contig or chromosome
        start: the 1-based start position of the Span (inclusive)
        end: the 1-based end position of the Span (inclusive)
        strand: the strand of the Span (POSITIVE by default)
    """

    refname: str
    start: int
    end: int
    strand: Strand = Strand.POSITIVE

    def __post_init__(self) -> None:
        if self.refname.strip() == "":
            raise ValueError("Ref name must be populated")
        if self.start < 1:
            raise ValueError(f"Start position must be >=1; received start={self.start}")
        if self.end < self.start:
            raise ValueError(f"End must be >= start; received start={self.start}, end={self.end}")

    @classmethod
    def from_string(cls, line: str) -> "Span":
        """Creates a Span from an input string.

        The input string should be delimited by a colon (":"). If strand information is missing
        after splitting the string, the strand is assumed to be positive.

        Args:
            line: input string

        Returns:
            Span: the Span object generated from the input string

        Raises:
            ValueError: if there are not at least 2 colon-delimited fields in string

        Example:
            >>> span_string = "chr1:1-10:+"
            >>> Span.from_string(span_string)
            Span(refname='chr1', start=1, end=10, strand=<Strand.POSITIVE: '+'>)
        """
        parts = line.strip().split(":")
        if len(parts) == 3:
            refname, range_, strand_symbol = parts
            try:
                strand = Strand(strand_symbol[0])
            except ValueError as e:
                raise ValueError(
                    "Did not find valid strand information; " f"received {strand_symbol}"
                ) from e
        elif len(parts) == 2:
            refname, range_ = parts
            strand = Strand.POSITIVE
        else:
            raise ValueError(
                f"Could not parse line (expected 2 or 3 colon-delimited fields), received {line}"
            )
        try:
            start, end = map(int, range_.split("-"))
        except ValueError as e:
            raise ValueError(
                f"Could not cast positional information to int; received {range_}"
            ) from e
        return Span(refname=refname, start=start, end=end, strand=strand)

    def __str__(self) -> str:
        return f"{self.refname}:{self.start}-{self.end}:{self.strand}"

    @property
    def length(self) -> int:
        """Get the length of the span/interval."""
        return self.end - self.start + 1

    def overlaps(self, other: "Span") -> bool:
        """Returns True if the spans overlap one another, False otherwise."""
        return (
            (self.refname == other.refname)
            and (self.start <= other.end)
            and (self.end >= other.start)
        )

    def length_of_overlap_with(self, other: "Span") -> int:
        """Returns the length of the region which overlaps the other span, or zero if there is
        no overlap"""
        overlap: int = 0
        if self.overlaps(other=other):
            start = max(self.start, other.start)
            end = min(self.end, other.end)
            overlap = end - start + 1
        return overlap

    def get_offset(self, position: int) -> int:
        """Returns a coordinate position that is relative to the current span.

        Args:
            position: the (genomic) position to provide relative coordinates for (1-based)

        Returns:
            A 0-based position relative to the Span

        Raises:
            ValueError: if the provided position is outside the coordinates of the Span

        Example:

        ```python
        >>> test_span = Span(refname="chr1", start=10, end=20, strand=Strand.POSITIVE)
        >>> print(test_span.get_offset(11))
        1

        ```
        """
        if position < self.start or position > self.end:
            raise ValueError(f"Position not within this span: {position}")
        return position - self.start

    def get_bedlike_coords(self) -> BedLikeCoords:
        """Returns the zero-based, open-ended coordinates (start and end) of the
         (1-based) Span, for using in bed-like formats.

         Returns:
             a BedLikeCoords instance containing the 0-based and open-ended coordinates of
             the Span

         Example:

        ```python
        >>> test_span = Span(refname="chr1", start=10, end=20, strand=Strand.POSITIVE)
        >>> print(test_span.get_bedlike_coords())
        BedLikeCoords(start=9, end=20)

        ```
        """
        return BedLikeCoords(self.start - 1, self.end)

    def get_subspan(
        self, offset: int, subspan_length: int, strand: Optional[Strand] = None
    ) -> "Span":
        """Returns a Span with absolute coords from an offset and a length.
        The strand of the new Span will be strand if given, otherwise it will be
        self.strand.

        Args:
            offset: the difference between the start position of the subspan and that
                of the current span (0-based)
            subspan_length: the length of the new span
            strand: the strand of the new span (if given)

        Returns:
            a new Span with objective coords derived from input offset and length

        Raises:
            ValueError: if the offset is less than 0
            ValueError: if the length of the subspan is 0 or less
            ValueError: if the offset is greater than the length of the source span

        Example:

        ```python
        >>> span = Span(refname="chr1", start=1000,end=2000, strand=Strand.POSITIVE)
        >>> test = span.get_subspan(offset=5, subspan_length=10)
        >>> print(test)
        chr1:1005-1014:+
        >>> print(test.length) #length = end - start + 1
        10

        ```
        """

        if offset < 0:
            raise ValueError(f"Offset must be > 0, received start={offset}")
        if subspan_length <= 0:
            raise ValueError(
                f"Length of a subspan must be positive, received length={subspan_length}"
            )
        if offset >= self.length:
            raise ValueError(
                "Offset of a relative subspan must be < source span length, "
                f"received offset={offset}, length={subspan_length}"
            )
        if offset + subspan_length > self.length:
            raise ValueError(
                f"End of sub-span must be within source span: source start={self.start}, "
                f"offset={offset}, sub-span length={subspan_length}"
            )

        absolute_start = self.start + offset
        absolute_end = coordmath.get_closed_end(absolute_start, subspan_length)
        strand = self.strand if strand is None else strand
        return replace(self, start=absolute_start, end=absolute_end, strand=strand)

    def contains(self, comparison: Self) -> bool:
        """Checks whether one span is wholly contained within another span.
        Returns `True` if one span contains the other, otherwise returns `False`.
        Does not use strand information (a span is considered to contain the other
        if they are adjacent in position but on opposite strands)."""
        return (
            self.refname == comparison.refname
            and self.start <= comparison.start
            and comparison.end <= self.end
        )

    def _to_tuple(self, seq_dict: SequenceDictionary) -> Tuple[int, int, int, int]:
        """Returns a tuple of reference index, start position, end position, and strand
        (0 forward, 1 reverse)"""
        ref_index = seq_dict.by_name(self.refname).index
        strand = 0 if self.strand == Strand.POSITIVE else 1
        return ref_index, self.start, self.end, strand

    @staticmethod
    def compare(this: "Span", that: "Span", seq_dict: SequenceDictionary) -> int:
        """Compares this span to that span, ordering references using the given sequence dictionary.

        Args:
            this: the first span
            that: the second span
            seq_dict: the sequence dictionary used to order references

        Returns:
            -1 if the first span is less than the second span, 0 if equal, 1 otherwise
        """
        left_tuple = this._to_tuple(seq_dict=seq_dict)
        right_tuple = that._to_tuple(seq_dict=seq_dict)
        if left_tuple < right_tuple:
            return -1
        elif left_tuple > right_tuple:
            return 1
        else:
            return 0
