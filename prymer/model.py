from dataclasses import dataclass
from dataclasses import replace
from enum import StrEnum
from enum import unique
from functools import cached_property
from typing import Any
from typing import Callable
from typing import Dict
from typing import Generic
from typing import Iterator
from typing import Optional
from typing import Self
from typing import Tuple
from typing import TypeVar

from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.sequence import gc_content
from fgpyo.sequence import longest_dinucleotide_run_length
from fgpyo.sequence import longest_homopolymer_length
from fgpyo.util.metric import Metric

Numeric = TypeVar("Numeric", int, float)


@dataclass(slots=True, frozen=True, init=True)
class MinOptMax(Generic[Numeric]):
    """Stores a minimum, an optimal, and a maximum value (either all ints or all floats).

    `min` must be less than `max`. `opt` should be greater than the `min`
     value and less than the `max` value.

    The three values can be either int or float values but all must be of the same type within one
    MinOptMax object (for example, `min` cannot be a float while `max` is an int).

    Examples of interacting with the `MinOptMax` class

    ```python
    >>> thresholds = MinOptMax(min=1.0, opt=2.0, max=4.0)
    >>> print(thresholds)
    (min:1.0, opt:2.0, max:4.0)
    >>> list(thresholds)
    [1.0, 2.0, 4.0]

    ```

    Attributes:
        min: the minimum value
        opt: the optimal value
        max: the maximum value

    Raises:
        ValueError: if min > max
        ValueError: if `min` is not less than `opt` and `opt` is not less than `max`
    """

    min: Numeric
    opt: Numeric
    max: Numeric

    def __post_init__(self) -> None:
        dtype = type(self.min)
        if not isinstance(self.max, dtype) or not isinstance(self.opt, dtype):
            raise TypeError(
                "Min, opt, and max must be the same type; "
                f"received min: {dtype}, opt: {type(self.opt)}, max: {type(self.max)}"
            )
        if self.min > self.max:
            raise ValueError(
                f"Min must be no greater than max; received min: {self.min}, max: {self.max}"
            )
        if not (self.min <= self.opt <= self.max):
            raise ValueError(
                "Arguments must satisfy min <= opt <= max: "
                f"received min: {self.min}, opt: {self.opt}, max: {self.max}"
            )

    def __iter__(self) -> Iterator[float]:
        """Returns an iterator of min, opt, and max"""
        return iter([self.min, self.opt, self.max])

    def __str__(self) -> str:
        """Returns a string representation of min, opt, and max"""
        return f"(min:{self.min}, opt:{self.opt}, max:{self.max})"


@unique
class Strand(StrEnum):
    """Represents the strand of a span to the genome."""

    POSITIVE = "+"
    NEGATIVE = "-"


@dataclass(init=True, frozen=True, eq=True)
class Span(Metric["Span"]):
    """Represents the span of some sequence (target, primer, amplicon, etc.) to a genome.

    Attributes:
        refname: name of a reference sequence or contig or chromosome
        start: the 1-based start position of the Span (inclusive)
        end: the 1-based end position of the Span (inclusive)
        strand: the strand of the Span (POSITIVE by default)

    Examples:

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

        ```python
        >>> span_string = "chr1:1-10:+"
        >>> Span.from_string(span_string)
        Span(refname='chr1', start=1, end=10, strand=<Strand.POSITIVE: '+'>)
        ```

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
        """Returns a string version of the span, e.g. `chr1:100-200:+`."""
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
        """Returns the relative offset of a 1-based genomic coordinate contained within the span.

        Args:
            position: the 1-based genomic position whose offset to calculate

        Returns:
            A 0-based offset relative to the Span

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
        absolute_end = absolute_start + subspan_length - 1
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


@dataclass(frozen=True, init=True, kw_only=True, slots=True)
class Oligo(Metric["Oligo"]):
    """Stores the properties of an oligo, which can be either a primer or internal oligo.

    The penalty score of the design is emitted by Primer3 and controlled by the corresponding
    design parameters. The penalty for a primer is set by the combination of
    `PrimerAndAmpliconParameters` and `PrimerWeights`. A probe penalty is set by `ProbeParameters`
    and `ProbeWeights`.

    The span of an `Oligo` object represents the mapping of the oligo to the genome. `Oligo` objects
    can optionally have a `tail` sequence to prepend to the 5' end of the primer or probe, as well
    as a `name` for downstream record keeping.

    `Oligo` objects can also store thermodynamic characteristics of primer and/or probe
    design. These include the calculated melting temperatures of the most stable hairpin structure,
    self-, and end- complementarity. These values can be emitted by Primer3.

    These thermodynamic characteristics are meant to quantify how likely it is the primer or probe
    will bind to itself (e.g., instead of the target DNA). A melting temperature close to or greater
    than the intended melting temperature for target binding indicates the primer or probe is
    likely to form stable hairpins or dimers, leading to reduced efficiency of the reaction.

    Examples of interacting with the `Oligo` class

    ```python
    >>> from prymer import Span, Strand
    >>> oligo_span = Span(refname="chr1", start=1, end=20)
    >>> oligo = Oligo(tm=70.0, penalty=-123.0, span=oligo_span, bases="AGCT" * 5)
    >>> oligo.longest_hp_length
    1
    >>> oligo.length
    20
    >>> oligo.name is None
    True
    >>> oligo = Oligo(tm=70.0, penalty=-123.0, span=oligo_span, bases="GACGG"*4)
    >>> oligo.longest_hp_length
    3
    >>> oligo.untailed_length
    20
    >>> oligo.tailed_length
    20
    >>> primer = oligo.with_tail(tail="GATTACA")
    >>> primer.untailed_length
    20
    >>> primer.tailed_length
    27
    >>> primer = primer.with_name(name="fwd_primer")
    >>> primer.name
    'fwd_primer'

    ```

    Oligos may also be written to a file and subsequently read back in, as the `Oligo` class is an
    `fgpyo` `Metric` class:

    ```python
    >>> from pathlib import Path
    >>> left_span = Span(refname="chr1", start=1, end=20)
    >>> left = Oligo(tm=70.0, penalty=-123.0, span=left_span, bases="G"*20)
    >>> right_span = Span(refname="chr1", start=101, end=120)
    >>> right = Oligo(tm=70.0, penalty=-123.0, span=right_span, bases="T"*20)
    >>> path = Path("/tmp/path/to/primers.txt")
    >>> Oligo.write(path, left, right)  # doctest: +SKIP
    >>> primers = Oligo.read(path)  # doctest: +SKIP
    >>> list(primers)  # doctest: +SKIP
    [
        Oligo(tm=70.0, penalty=-123.0, span=amplicon_span, bases="G"*20),
        Oligo(tm=70.0, penalty=-123.0, span=amplicon_span, bases="T"*20)
    ]

    ```

    Attributes:
        bases: the base sequence of the oligo (excluding any tail)
        tm: the calculated melting temperature of the oligo
        span: the mapping of the primer to the genome
        penalty: the penalty or score for the oligo
        name: an optional name to use for the primer
        tm_homodimer: calculated melting temperature that represents
            the tendency of an oligo to bind to itself (self-complementarity)
        tm_3p_anchored_homodimer: calculated melting temperature that represents
            the tendency of a primer to bind to the 3'-END of an identical primer
        tm_secondary_structure: calculated melting temperature of the oligo hairpin structure
        tail: an optional tail sequence to put on the 5' end of the primer

    """

    bases: str
    tm: float
    span: Span
    penalty: float
    name: Optional[str] = None
    tm_homodimer: Optional[float] = None
    tm_3p_anchored_homodimer: Optional[float] = None
    tm_secondary_structure: Optional[float] = None
    tail: Optional[str] = None

    def __post_init__(self) -> None:
        if len(self.bases) == 0:
            raise ValueError("Bases must not be an empty string.")
        elif self.span.length != len(self.bases):
            raise ValueError(
                f"Conflicting lengths: span={self.span.length}, bases={len(self.bases)}"
            )

    @property
    def id(self) -> str:
        """Returns the name if there is one, otherwise the span."""
        return self.name if self.name is not None else f"{self.span}"

    @property
    def length(self) -> int:
        """Length of un-tailed oligo."""
        return self.span.length

    @property
    def untailed_length(self) -> int:
        """Length of un-tailed oligo."""
        return self.span.length

    @property
    def tailed_length(self) -> int:
        """Length of tailed oligo."""
        return self.span.length + (0 if self.tail is None else len(self.tail))

    @cached_property
    def longest_hp_length(self) -> int:
        """Length of longest homopolymer in the oligo."""
        return longest_homopolymer_length(self.bases)

    @cached_property
    def percent_gc_content(self) -> float:
        """The GC of the amplicon sequence in the range 0-100."""
        return round(gc_content(self.bases) * 100, 3)

    @cached_property
    def longest_dinucleotide_run_length(self) -> int:
        """
        Number of bases in the longest dinucleotide run in a oligo.

        A dinucleotide run is when length two repeat-unit is repeated. For example,
        TCTC (length = 4) or ACACACACAC (length = 10). If there are no such runs, returns 2
        (or 0 if there are fewer than 2 bases).
        """
        return longest_dinucleotide_run_length(self.bases)

    def with_tail(self, tail: str) -> "Oligo":
        """Returns a copy of the oligo with the tail sequence attached."""
        return replace(self, tail=tail)

    def with_name(self, name: str) -> "Oligo":
        """Returns a copy of oligo object with the given name."""
        return replace(self, name=name)

    @property
    def bases_with_tail(self) -> str:
        """
        Returns the sequence of the oligo prepended by the tail.

        If `tail` is None, only return `bases`.
        """
        if self.tail is None:
            return self.bases
        else:
            return f"{self.tail}{self.bases}"

    def __str__(self) -> str:
        """Returns a string representation of this oligo."""
        return f"{self.bases}\t{self.tm:.2f}\t{self.penalty:.2f}\t{self.span}"

    @classmethod
    def _parsers(cls) -> Dict[type, Callable[[str], Any]]:
        return {
            Span: lambda value: Span.from_string(value),
        }

    @staticmethod
    def compare(this: "Oligo", that: "Oligo", seq_dict: SequenceDictionary) -> int:
        """Compares this oligo to that oligo by their span, ordering references using the given
        sequence dictionary.

        Args:
            this: the first oligo
            that: the second oligo
            seq_dict: the sequence dictionary used to order references

        Returns:
            -1 if this oligo is less than the that oligo, 0 if equal, 1 otherwise
        """
        return Span.compare(this=this.span, that=that.span, seq_dict=seq_dict)


@dataclass(frozen=True, init=True, kw_only=True)
class PrimerPair:
    """
    Represents a pair of primers that work together to amplify an amplicon. The coordinates of the
    amplicon are determined by the start of the left primer and the end of the right primer.

    Optional attributes include naming information to keep track of the pair and design parameters
    used to create the pairing.

    Attributes:
        left_primer: the left primer in the pair
        right_primer: the right primer in the pair
        amplicon_sequence: an optional sequence of the expected amplicon
        amplicon_tm: the melting temperature of the expected amplicon
        penalty: the penalty value assigned by primer3 to the primer pair
        name: an optional name to use for the primer pair

    Raises:
        ValueError: if the chromosomes of the left and right primers are not the same

    Examples:
        ```python
        >>> from prymer import Strand
        >>> left_span = Span(refname="chr1", start=1, end=20)
        >>> left_primer = Oligo(tm=70.0, penalty=-123.0, span=left_span, bases="G"*20)
        >>> right_span = Span(refname="chr1", start=101, end=120, strand=Strand.NEGATIVE)
        >>> right_primer = Oligo(tm=70.0, penalty=-123.0, span=right_span, bases="T"*20)
        >>> primer_pair = PrimerPair( \
            left_primer=left_primer, \
            right_primer=right_primer, \
            amplicon_sequence=None, \
            amplicon_tm=70.0, \
            penalty=-123.0, \
            name="foobar" \
        )
        >>> primer_pair.amplicon
        Span(refname='chr1', start=1, end=120, strand=<Strand.POSITIVE: '+'>)
        >>> primer_pair.span
        Span(refname='chr1', start=1, end=120, strand=<Strand.POSITIVE: '+'>)
        >>> list(primer_pair)
        [Oligo(bases='GGGGGGGGGGGGGGGGGGGG', tm=70.0, span=Span(refname='chr1', start=1, end=20, strand=<Strand.POSITIVE: '+'>), penalty=-123.0, name=None, tm_homodimer=None, tm_3p_anchored_homodimer=None, tm_secondary_structure=None, tail=None), Oligo(bases='TTTTTTTTTTTTTTTTTTTT', tm=70.0, span=Span(refname='chr1', start=101, end=120, strand=<Strand.NEGATIVE: '-'>), penalty=-123.0, name=None, tm_homodimer=None, tm_3p_anchored_homodimer=None, tm_secondary_structure=None, tail=None)]
    ```

    """  # noqa: E501

    left_primer: Oligo
    right_primer: Oligo
    amplicon_tm: float
    penalty: float
    name: Optional[str] = None
    amplicon_sequence: Optional[str] = None

    def __post_init__(self) -> None:
        # Force generation of the amplicon span, which validates the relative positions
        # of the left and right primers
        amp = self.amplicon

        # If supplied, bases must be the same length as the span length
        if self.bases is not None:
            if len(self.bases) == 0:
                raise ValueError("Bases must not be an empty string")
            elif len(self.bases) != amp.length:
                raise ValueError(
                    f"Conflicting lengths: span={amp.length}, sequence={len(self.bases)}"
                )

    @property
    def id(self) -> str:
        """Returns the name if there is one, otherwise the span."""
        return self.name if self.name is not None else str(self.span)

    @cached_property
    def amplicon(self) -> Span:
        """Returns the mapping for the amplicon"""
        return self.calculate_amplicon_span(self.left_primer, self.right_primer)

    @property
    def span(self) -> Span:
        """Returns the mapping for the amplicon"""
        return self.amplicon

    @property
    def bases(self) -> Optional[str]:
        """Returns the bases of the amplicon sequence"""
        return self.amplicon_sequence

    @property
    def length(self) -> int:
        """Returns the length of the amplicon"""
        return self.amplicon.length

    @cached_property
    def percent_gc_content(self) -> float:
        """The GC of the amplicon sequence in the range 0-100, or 0 if amplicon sequence is None."""
        if self.bases is None:
            return 0.0
        else:
            return round(gc_content(self.bases) * 100, 3)

    def with_tails(self, left_tail: str, right_tail: str) -> "PrimerPair":
        """
        Returns a copy of the primer pair where the left and right primers are tailed.

        Args:
            left_tail: The tail to add to the left primer
            right_tail: The tail to add to the right primer

        Returns:
            A copy of the primer pair with the tail(s) added to the primers
        """
        return replace(
            self,
            left_primer=self.left_primer.with_tail(left_tail),
            right_primer=self.right_primer.with_tail(right_tail),
        )

    def with_names(self, pp_name: str, lp_name: str, rp_name: str) -> "PrimerPair":
        """
        Returns a copy of the primer pair with names assigned to the primer pair,
        left primer, and right primer.

        Args:
            pp_name: The optional name of the primer pair
            lp_name: The optional name of the left primer
            rp_name: The optional name of the right primer

        Returns:
            A copy of the primer pair with the provided names assigned
        """
        return replace(
            self,
            name=pp_name,
            left_primer=self.left_primer.with_name(lp_name),
            right_primer=self.right_primer.with_name(rp_name),
        )

    def __iter__(self) -> Iterator[Oligo]:
        """Returns an iterator of left and right primers"""
        return iter([self.left_primer, self.right_primer])

    def __str__(self) -> str:
        """Returns a string representation of the primer pair"""
        sequence = self.amplicon_sequence if self.amplicon_sequence else "*"
        return (
            f"{self.left_primer}\t{self.right_primer}\t{sequence}\t"
            + f"{self.amplicon_tm}\t{self.penalty}"
        )

    @staticmethod
    def calculate_amplicon_span(left_primer: Oligo, right_primer: Oligo) -> Span:
        """
        Calculates the amplicon Span from the left and right primers.

        Args:
            left_primer: the left primer for the amplicon
            right_primer: the right primer for the amplicon

        Returns:
            a Span starting at the first base of the left primer and ending at the last base of
            the right primer

        Raises:
            ValueError: If `left_primer` and `right_primer` have different reference names.
            ValueError: If `left_primer` doesn't start before the right primer.
            ValueError: If `right_primer` ends before `left_primer`.
        """
        # Require that `left_primer` and `right_primer` both map to the same reference sequence
        if left_primer.span.refname != right_primer.span.refname:
            raise ValueError(
                "Left and right primers are on different references. "
                f"Left primer ref: {left_primer.span.refname}. "
                f"Right primer ref: {right_primer.span.refname}"
            )

        # Require that the left primer starts before the right primer
        if left_primer.span.start > right_primer.span.start:
            raise ValueError(
                "Left primer does not start before the right primer. "
                f"Left primer span: {left_primer.span}, "
                f"Right primer span: {right_primer.span}"
            )

        # Require that the left primer starts before the right primer
        if right_primer.span.end < left_primer.span.end:
            raise ValueError(
                "Right primer ends before left primer ends. "
                f"Left primer span: {left_primer.span}, "
                f"Right primer span: {right_primer.span}"
            )

        return Span(left_primer.span.refname, left_primer.span.start, right_primer.span.end)

    @staticmethod
    def compare(
        this: "PrimerPair",
        that: "PrimerPair",
        seq_dict: SequenceDictionary,
        by_amplicon: bool = True,
    ) -> int:
        """Compares this primer pair to that primer pair by their span, ordering references using
        the given sequence dictionary.

        Args:
            this: the first primer pair
            that: the second primer pair
            seq_dict: the sequence dictionary used to order references
            by_amplicon: ture to compare using the amplicon property on a primer pair, false to
                compare first using the left primer then the right primer

        Returns:
            -1 if this primer pair is less than the that primer pair, 0 if equal, 1 otherwise
        """
        if by_amplicon:
            return Span.compare(this=this.amplicon, that=that.amplicon, seq_dict=seq_dict)
        else:
            retval = Oligo.compare(this=this.left_primer, that=that.left_primer, seq_dict=seq_dict)
            if retval == 0:
                retval = Oligo.compare(
                    this=this.right_primer, that=that.right_primer, seq_dict=seq_dict
                )
            return retval
