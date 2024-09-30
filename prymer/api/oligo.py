"""
# Oligo Class and Methods

This module contains a class and class methods to represent an oligo designed by Primer3.

Oligos can represent single primer and/or internal probe designs.

Class attributes include the base sequence, melting temperature, and the score of the oligo. The
mapping of the oligo to the genome is also stored.

Optional attributes include naming information and a tail sequence to attach to the 5' end of the
oligo (if applicable). Optional attributes also include the thermodynamic results from Primer3.

## Examples of interacting with the `Oligo` class

```python
>>> from prymer.api.span import Span, Strand
>>> oligo_span = Span(refname="chr1", start=1, end=20)
>>> oligo = Oligo(tm=70.0, penalty=-123.0, span=oligo_span, bases="AGCT" * 5)
>>> oligo.longest_hp_length()
1
>>> oligo.length
20
>>> oligo.name is None
True
>>> oligo = Oligo(tm=70.0, penalty=-123.0, span=oligo_span, bases="GACGG"*4)
>>> oligo.longest_hp_length()
3
>>> oligo.untailed_length()
20
>>> oligo.tailed_length()
20
>>> primer = oligo.with_tail(tail="GATTACA")
>>> primer.untailed_length()
20
>>> primer.tailed_length()
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

from prymer.api.oligo_like import MISSING_BASES_STRING
from prymer.api.oligo_like import OligoLike
from prymer.api.span import Span


@dataclass(frozen=True, init=True, kw_only=True, slots=True)
class Oligo(OligoLike, Metric["Oligo"]):
    """Stores the properties of the designed oligo.

    Oligos can include both single primer and internal probe designs. The penalty score of the
    design is emitted by Primer3 and controlled by the corresponding design parameters.
    The penalty for a primer is set by the combination of `PrimerAndAmpliconParameters` and
    `PrimerWeights`, whereas a probe penalty is set by `ProbeParameters` and `ProbeWeights`.

    The values for `self_any`, `self_any_th`, `self_end`, `self_end_th`, and `hairpin_th`
    are emitted by Primer3 as part of oligo design. These attributes are optional to maintain
    flexibility for reading in and writing `Oligo` objects, espeically when design settings are
    inconsistent.

    Attributes:
        tm: the calculated melting temperature of the oligo
        penalty: the penalty or score for the oligo
        span: the mapping of the primer to the genome
        self_any: probe self-complementarity, expressed as local alignment score
        self_any_th: probe self-complementarity, expressed as melting temperature
        self_end: 3' end complementarity, expressed as local alignment score
        self_end_th: 3' end complementarity, expressed as melting temperature
        hairpin_th: hairpin formation thermodynamics of the oligo as calculated by Primer3
        bases: the base sequence of the oligo (excluding any tail)
        tail: an optional tail sequence to put on the 5' end of the primer
        name: an optional name to use for the primer

    """

    tm: float
    penalty: float
    span: Span
    self_any: Optional[float] = None
    self_any_th: Optional[float] = None
    self_end: Optional[float] = None
    self_end_th: Optional[float] = None
    hairpin_th: Optional[float] = None
    bases: Optional[str] = None
    tail: Optional[str] = None

    def __post_init__(self) -> None:
        super(Oligo, self).__post_init__()

    def longest_hp_length(self) -> int:
        """Length of longest homopolymer in the oligo."""
        if self.bases is None:
            return 0
        else:
            return longest_homopolymer_length(self.bases)

    @property
    def length(self) -> int:
        """Length of un-tailed oligo."""
        return self.span.length

    def untailed_length(self) -> int:
        """Length of un-tailed oligo."""
        return self.span.length

    def tailed_length(self) -> int:
        """Length of tailed oligo."""
        return self.span.length if self.tail is None else self.span.length + len(self.tail)

    def longest_dinucleotide_run_length(self) -> int:
        """Number of bases in the longest dinucleotide run in a oligo.

        A dinucleotide run is when length two repeat-unit is repeated. For example,
        TCTC (length = 4) or ACACACACAC (length = 10). If there are no such runs, returns 2
        (or 0 if there are fewer than 2 bases)."""
        return longest_dinucleotide_run_length(self.bases)

    def with_tail(self, tail: str) -> "Oligo":
        """Returns a copy of the oligo with the tail sequence attached."""
        return replace(self, tail=tail)

    def with_name(self, name: str) -> "Oligo":
        """Returns a copy of oligo object with the given name."""
        return replace(self, name=name)

    def bases_with_tail(self) -> Optional[str]:
        """
        Returns the sequence of the oligo prepended by the tail.

        If `tail` is None, only return `bases`.
        """
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
        Returns a string representation of this oligo
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
