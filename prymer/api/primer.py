"""
# Primer Class and Methods

This module contains a class and class methods to represent a primer (e.g. designed by Primer3)

Class attributes include the primer sequence, melting temperature, and the score of the primer. The
mapping of the primer to the genome is also stored.

Optional attributes include naming information and a tail sequence  to attach to the 5' end of the
primer (if applicable).

## Examples of interacting with the `Primer` class

```python
>>> from prymer.api.span import Span, Strand
>>> primer_span = Span(refname="chr1", start=1, end=20)
>>> primer = Primer(tm=70.0, penalty=-123.0, span=primer_span)
>>> primer.longest_hp_length()
0
>>> primer.length
20
>>> primer.name is None
True
>>> primer = Primer(tm=70.0, penalty=-123.0, span=primer_span, bases="GACGG"*4)
>>> primer.longest_hp_length()
3
>>> primer.untailed_length()
20
>>> primer.tailed_length()
20
>>> primer = primer.with_tail(tail="GATTACA")
>>> primer.untailed_length()
20
>>> primer.tailed_length()
27
>>> primer = primer.with_name(name="foobar")
>>> primer.name
'foobar'

```

Primers may also be written to a file and subsequently read back in, as the `Primer` class is an
`fgpyo` `Metric` class:

```python
>>> from pathlib import Path
>>> left_span = Span(refname="chr1", start=1, end=20)
>>> left = Primer(tm=70.0, penalty=-123.0, span=left_span, bases="G"*20)
>>> right_span = Span(refname="chr1", start=101, end=120)
>>> right = Primer(tm=70.0, penalty=-123.0, span=right_span, bases="T"*20)
>>> path = Path("/tmp/path/to/primers.txt")
>>> Primer.write(path, left, right)  # doctest: +SKIP
>>> primers = Primer.read(path)  # doctest: +SKIP
>>> list(primers)  # doctest: +SKIP
[
    Primer(tm=70.0, penalty=-123.0, span=amplicon_span, bases="G"*20),
    Primer(tm=70.0, penalty=-123.0, span=amplicon_span, bases="T"*20)
]

```
"""

from dataclasses import dataclass
from dataclasses import replace
from typing import Any
from typing import Callable
from typing import Dict
from typing import Optional

from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.sequence import longest_dinucleotide_run_length
from fgpyo.sequence import longest_homopolymer_length
from fgpyo.util.metric import Metric

from prymer.api.primer_like import MISSING_BASES_STRING
from prymer.api.primer_like import PrimerLike
from prymer.api.span import Span


@dataclass(frozen=True, init=True, kw_only=True, slots=True)
class Primer(PrimerLike, Metric["Primer"]):
    """Stores the properties of the designed Primer.

    Attributes:
        bases: the base sequence of the primer (excluding any tail)
        tm: the calculated melting temperature of the primer
        penalty: the penalty or score for the primer
        span: the mapping of the primer to the genome
        name: an optional name to use for the primer
        tail: an optional tail sequence to put on the 5' end of the primer

    Example:
        #TODO
        <.....>
    """

    tm: float
    penalty: float
    span: Span
    bases: Optional[str] = None
    tail: Optional[str] = None

    def __post_init__(self) -> None:
        super(Primer, self).__post_init__()

    def longest_hp_length(self) -> int:
        """Length of longest homopolymer in the primer."""
        if self.bases is None:
            return 0
        else:
            return longest_homopolymer_length(self.bases)

    @property
    def length(self) -> int:
        """Length of un-tailed primer."""
        return self.span.length

    def untailed_length(self) -> int:
        """Length of un-tailed primer."""
        return self.span.length

    def tailed_length(self) -> int:
        """Length of tailed primer."""
        return self.span.length if self.tail is None else self.span.length + len(self.tail)

    def longest_dinucleotide_run_length(self) -> int:
        """Number of bases in the longest dinucleotide run in a primer.

        A dinucleotide run is when length two repeat-unit is repeated. For example,
        TCTC (length = 4) or ACACACACAC (length = 10). If there are no such runs, returns 2
        (or 0 if there are fewer than 2 bases)."""
        return longest_dinucleotide_run_length(self.bases)

    def with_tail(self, tail: str) -> "Primer":
        """Returns a copy of the primer with the tail sequence attached."""
        return replace(self, tail=tail)

    def with_name(self, name: str) -> "Primer":
        """Returns copy of primer object with the given name."""
        return replace(self, name=name)

    def bases_with_tail(self) -> Optional[str]:
        """
        Returns the sequence of the primer prepended by the tail.

        If either `bases` or `tail` are None, they shall be excluded. Return None if both are None.
        """
        if self.bases is None:
            return None if self.tail is None else self.tail
        if self.tail is None:
            return self.bases
        return f"{self.tail}{self.bases}"

    def to_bed12_row(self) -> str:
        """Returns the BED detail format view:
        https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7"""
        bed_coord = self.span.get_bedlike_coords()
        return "\t".join(
            map(
                str,
                [
                    self.span.refname,  # contig
                    bed_coord.start,  # start
                    bed_coord.end,  # end
                    self.id,  # name
                    500,  # score
                    self.span.strand.value,  # strand
                    bed_coord.start,  # thick start
                    bed_coord.end,  # thick end
                    "100,100,100",  # color
                    1,  # block count
                    f"{self.length}",  # block sizes
                    "0",  # block starts (relative to `start`)
                ],
            )
        )

    def __str__(self) -> str:
        """
        Returns a string representation of this primer
        """
        # If the bases field is None, replace with MISSING_BASES_STRING
        bases: str = self.bases if self.bases is not None else MISSING_BASES_STRING
        return f"{bases}\t{self.tm}\t{self.penalty}\t{self.span}"

    @classmethod
    def _parsers(cls) -> Dict[type, Callable[[str], Any]]:
        return {
            Span: lambda value: Span.from_string(value),
        }

    @staticmethod
    def compare(this: "Primer", that: "Primer", seq_dict: SequenceDictionary) -> int:
        """Compares this primer to that primer by their span, ordering references using the given
        sequence dictionary.

        Args:
            this: the first primer
            that: the second primer
            seq_dict: the sequence dictionary used to order references

        Returns:
            -1 if this primer is less than the that primer, 0 if equal, 1 otherwise
        """
        return Span.compare(this=this.span, that=that.span, seq_dict=seq_dict)
