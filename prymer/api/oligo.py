"""
# Oligo Class and Methods

This module contains a class and class methods to represent an oligo (e.g., designed by Primer3).

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
"""

from dataclasses import dataclass
from dataclasses import replace
from functools import cached_property
from typing import Any
from typing import Callable
from typing import Dict
from typing import Optional

from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.sequence import gc_content
from fgpyo.sequence import longest_dinucleotide_run_length
from fgpyo.sequence import longest_homopolymer_length
from fgpyo.util.metric import Metric

from prymer.api.span import Span


@dataclass(frozen=True, init=True, kw_only=True, slots=True)
class Oligo(Metric["Oligo"]):
    """Stores the properties of the designed oligo.

    Oligos can include both single primer and internal probe designs. Primer3 emits information
    about the specific properties of a primer and/or probe, including the base sequence,
    melting temperature (tm), and score.

    The penalty score of the design is emitted by Primer3 and controlled by the corresponding
    design parameters. The penalty for a primer is set by the combination of
    `PrimerAndAmpliconParameters` and `PrimerWeights`.
    A probe penalty is set by `ProbeParameters` and `ProbeWeights`.

    The span of an `Oligo` object represents the mapping of the oligo
    to the genome. `Oligo` objects can optionally have a `tail` sequence to prepend to the 5' end
    of the primer or probe, as well as a `name` for downstream record keeping.

    `Oligo` objects can also store thermodynamic characteristics of primer and/or probe
    design. These include the calculated melting temperatures of the most stable hairpin structure,
    self-, and end- complementarity. These values can be emitted by Primer3.

    These thermodynamic characteristics are meant to quantify how likely it is the primer or probe
    will bind to itself (e.g., instead of the target DNA). A melting temperature close to or greater
    than the intended melting temperature for target binding indicates the primer or probe is
    likely to form stable hairpins or dimers, leading to reduced efficiency of the reaction.

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
