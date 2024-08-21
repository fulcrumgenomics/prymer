"""
# Primer Pair Class and Methods

This module contains the [`PrimerPair`][prymer.api.primer_pair.PrimerPair] class and
class methods to represent a primer pair.  The primer pair is comprised of a left and right primer
that work together to amplify an amplicon.

Class attributes include each of the primers (represented by a
[`Primer`][prymer.api.primer.Primer] object), information about the expected amplicon
(positional information about how the amplicon maps to the genome, the sequence, and its melting
temperature), as well as a score of the primer pair (e.g. as emitted by Primer3).

Optional attributes include naming information to keep track of the pair and design parameters
used to create the pairing.

## Examples

```python
>>> from prymer.api.span import Strand
>>> left_span = Span(refname="chr1", start=1, end=20)
>>> left_primer = Primer(tm=70.0, penalty=-123.0, span=left_span, bases="G"*20)
>>> right_span = Span(refname="chr1", start=101, end=120, strand=Strand.NEGATIVE)
>>> right_primer = Primer(tm=70.0, penalty=-123.0, span=right_span, bases="T"*20)
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
>>> primer_pair.inner
Span(refname='chr1', start=21, end=100, strand=<Strand.POSITIVE: '+'>)

>>> list(primer_pair)
[Primer(name=None, tm=70.0, penalty=-123.0, span=Span(refname='chr1', start=1, end=20, strand=<Strand.POSITIVE: '+'>), bases='GGGGGGGGGGGGGGGGGGGG', tail=None), Primer(name=None, tm=70.0, penalty=-123.0, span=Span(refname='chr1', start=101, end=120, strand=<Strand.NEGATIVE: '-'>), bases='TTTTTTTTTTTTTTTTTTTT', tail=None)]

```
"""  # noqa: E501

from dataclasses import dataclass
from dataclasses import field
from dataclasses import replace
from typing import Iterator
from typing import Optional

from fgpyo.fasta.sequence_dictionary import SequenceDictionary

from prymer.api.primer import Primer
from prymer.api.primer_like import MISSING_BASES_STRING
from prymer.api.primer_like import PrimerLike
from prymer.api.span import Span


@dataclass(frozen=True, init=True, kw_only=True, slots=True)
class PrimerPair(PrimerLike):
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

    left_primer: Primer
    right_primer: Primer
    amplicon_tm: float
    penalty: float
    amplicon_sequence: Optional[str] = None
    _amplicon: Span = field(init=False)

    def __post_init__(self) -> None:
        # Derive the amplicon from the left and right primers. This must be done before
        # calling super() as `PrimerLike.id` depends on the amplicon being set
        object.__setattr__(self, "_amplicon", self._calculate_amplicon())
        super(PrimerPair, self).__post_init__()

    @property
    def amplicon(self) -> Span:
        """Returns the mapping for the amplicon"""
        return self._amplicon

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

    @property
    def inner(self) -> Span:
        """
        Returns the inner region of the amplicon (not including the primers). I.e. the region of
        the genome covered by the primer pair, without the primer regions.  If the primers
        overlap, then the inner mapping is the midpoint at where they overlap
        """
        if self.left_primer.span.overlaps(self.right_primer.span):
            # Use a flooring division as these values are all ints
            midpoint = (self.left_primer.span.end + self.right_primer.span.start) // 2
            return replace(
                self.left_primer.span,
                start=midpoint,
                end=midpoint,
            )
        else:
            return replace(
                self.left_primer.span,
                start=self.left_primer.span.end + 1,
                end=self.right_primer.span.start - 1,
            )

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

    def to_bed12_row(self) -> str:
        """
        Returns the BED detail format view: https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7.

        NB: BED is 0-based and Prymer is 1-based, so we need to convert
        """
        block_sizes = ",".join(
            [
                f"{gm.length}"
                for gm in [
                    self.left_primer.span,
                    self.inner,
                    self.right_primer.span,
                ]
            ]
        )

        block_starts = ",".join(
            [
                f"{self.amplicon.get_offset(gm.start)}"
                for gm in [
                    self.left_primer.span,
                    self.inner,
                    self.right_primer.span,
                ]
            ],
        )
        bed_like_coords = self.span.get_bedlike_coords()
        return "\t".join(
            map(
                str,
                [
                    self.span.refname,  # contig
                    bed_like_coords.start,  # start
                    bed_like_coords.end,  # end
                    self.id,  # name
                    500,  # score
                    self.span.strand.value,  # strand
                    bed_like_coords.start,  # thick start
                    bed_like_coords.end,  # thick end
                    "100,100,100",  # color
                    3,  # block count
                    block_sizes,
                    block_starts,  # relative to `start`
                ],
            )
        )

    def __iter__(self) -> Iterator[Primer]:
        """Returns an iterator of left and right primers"""
        return iter([self.left_primer, self.right_primer])

    def __str__(self) -> str:
        """Returns a string representation of the primer pair"""
        sequence = self.amplicon_sequence if self.amplicon_sequence else MISSING_BASES_STRING
        return (
            f"{self.left_primer}\t{self.right_primer}\t{sequence}\t"
            + f"{self.amplicon_tm}\t{self.penalty}"
        )

    def _calculate_amplicon(self) -> Span:
        """
        Calculates the amplicon from the left and right primers, spanning from the start of the
        left primer to the end of the right primer.
        """

        # Require that `left_primer` and `right_primer` both map to the same reference sequence
        if self.left_primer.span.refname != self.right_primer.span.refname:
            raise ValueError(
                "The reference must be the same across primers in a pair; received "
                f"left primer ref: {self.left_primer.span.refname}, "
                f"right primer ref: {self.right_primer.span.refname}"
            )

        # Require that the left primer does not start to the right of the right primer
        if self.left_primer.span.start > self.right_primer.span.end:
            raise ValueError(
                "Left primer start must be less than or equal to right primer end; received "
                "left primer genome span: {self.left_primer.span}, "
                "right primer genome span: {self.right_primer.span}"
            )

        return replace(self.left_primer.span, end=self.right_primer.span.end)

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
            retval = Primer.compare(this=this.left_primer, that=that.left_primer, seq_dict=seq_dict)
            if retval == 0:
                retval = Primer.compare(
                    this=this.right_primer, that=that.right_primer, seq_dict=seq_dict
                )
            return retval
