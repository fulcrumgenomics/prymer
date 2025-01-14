"""
# Primer Pair Class and Methods

This module contains the [`PrimerPair`][prymer.api.primer_pair.PrimerPair] class and
class methods to represent a primer pair.  The primer pair is comprised of a left and right primer
that work together to amplify an amplicon.

Class attributes include each of the primers (represented by an
[`Oligo`][prymer.api.primer.Oligo] object), information about the expected amplicon
(positional information about how the amplicon maps to the genome, the sequence, and its melting
temperature), as well as a score of the primer pair (e.g. as emitted by Primer3).

Optional attributes include naming information to keep track of the pair and design parameters
used to create the pairing.

## Examples

```python
>>> from prymer.api.span import Strand
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
"""  # noqa: E501

from dataclasses import dataclass
from dataclasses import replace
from functools import cached_property
from typing import Iterator
from typing import Optional

from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.sequence import gc_content

from prymer.api.oligo import Oligo
from prymer.api.span import Span


@dataclass(frozen=True, init=True, kw_only=True)
class PrimerPair:
    """
    Represents a pair of primers that work together to amplify an amplicon. The
    coordinates of the amplicon are determined to span from the start of the left
    primer through the end of the right primer.

    Attributes:
        left_primer: the left primer in the pair
        right_primer: the right primer in the pair
        amplicon_sequence: an optional sequence of the expected amplicon
        amplicon_tm: the melting temperature of the expected amplicon
        penalty: the penalty value assigned by primer3 to the primer pair
        name: an optional name to use for the primer pair

    Raises:
        ValueError: if the chromosomes of the left and right primers are not the same
    """

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
