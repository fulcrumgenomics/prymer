"""
# Primer3Parameters Class and Methods

The [`Primer3Parameters`][prymer.primer3.primer3_parameters.Primer3Parameters] class stores
user input and maps it to the correct Primer3 fields.

Primer3 considers many criteria for primer design, including characteristics of candidate primers
and the resultant amplicon product, as well as potential complications (off-target priming,
primer dimer formation). Users can specify many of these constraints in Primer3,
some of which are used to quantify a "score" for each primer design.

The Primer3Parameters class stores commonly used constraints for primer design: GC content, melting
temperature, and size of both primers and expected amplicon. Additional criteria include the maximum
homopolymer length, ambiguous bases, and bases in a dinucleotide run within a primer. By default,
primer design avoids masked bases, returns 5 primers, and sets the GC clamp to be no larger than 5.

The `to_input_tags()` method in `Primer3Parameters` converts these parameters into tag-values pairs
for use when executing `Primer3`.

## Examples

```python
>>> params = Primer3Parameters( \
    amplicon_sizes=MinOptMax(min=100, max=250, opt=200), \
    amplicon_tms=MinOptMax(min=55.0, max=100.0, opt=70.0), \
    primer_sizes=MinOptMax(min=29, max=31, opt=30), \
    primer_tms=MinOptMax(min=63.0, max=67.0, opt=65.0), \
    primer_gcs=MinOptMax(min=30.0, max=65.0, opt=45.0), \
)
>>> for tag, value in params.to_input_tags().items(): \
    print(f"{tag.value} -> {value}")
PRIMER_NUM_RETURN -> 5
PRIMER_PRODUCT_OPT_SIZE -> 200
PRIMER_PRODUCT_SIZE_RANGE -> 100-250
PRIMER_PRODUCT_MIN_TM -> 55.0
PRIMER_PRODUCT_OPT_TM -> 70.0
PRIMER_PRODUCT_MAX_TM -> 100.0
PRIMER_MIN_SIZE -> 29
PRIMER_OPT_SIZE -> 30
PRIMER_MAX_SIZE -> 31
PRIMER_MIN_TM -> 63.0
PRIMER_OPT_TM -> 65.0
PRIMER_MAX_TM -> 67.0
PRIMER_MIN_GC -> 30.0
PRIMER_OPT_GC_PERCENT -> 45.0
PRIMER_MAX_GC -> 65.0
PRIMER_GC_CLAMP -> 0
PRIMER_MAX_END_GC -> 5
PRIMER_MAX_POLY_X -> 5
PRIMER_MAX_NS_ACCEPTED -> 1
PRIMER_LOWERCASE_MASKING -> 1

```
"""

from dataclasses import dataclass
from typing import Any

from prymer.api.minoptmax import MinOptMax
from prymer.primer3.primer3_input_tag import Primer3InputTag


@dataclass(frozen=True, init=True, slots=True)
class Primer3Parameters:
    """Holds common primer design options that Primer3 uses to inform primer design.

    Attributes:
        amplicon_sizes: the min, optimal, and max amplicon size
        amplicon_tms: the min, optimal, and max amplicon melting temperatures
        primer_sizes: the min, optimal, and max primer size
        primer_tms: the min, optimal, and max primer melting temperatures
        primer_gcs: the min and maximal GC content for individual primers
        gc_clamp: the min and max number of Gs and Cs in the 3' most N bases
        primer_max_polyX: the max homopolymer length acceptable within a primer
        primer_max_Ns: the max number of ambiguous bases acceptable within a primer
        primer_max_dinuc_bases: the maximal number of bases in a dinucleotide run in a primer
        avoid_masked_bases: whether Primer3 should avoid designing primers in soft-masked regions
        number_primers_return: the number of primers to return

    Please see the Primer3 manual for additional details: https://primer3.org/manual.html#globalTags

    """

    amplicon_sizes: MinOptMax[int]
    amplicon_tms: MinOptMax[float]
    primer_sizes: MinOptMax[int]
    primer_tms: MinOptMax[float]
    primer_gcs: MinOptMax[float]
    gc_clamp: tuple[int, int] = (0, 5)
    primer_max_polyX: int = 5
    primer_max_Ns: int = 1
    primer_max_dinuc_bases: int = 6
    avoid_masked_bases: bool = True
    number_primers_return: int = 5

    def __post_init__(self) -> None:
        if self.primer_max_dinuc_bases % 2 == 1:
            raise ValueError("Primer Max Dinuc Bases must be an even number of bases")
        if not isinstance(self.amplicon_sizes.min, int) or not isinstance(
            self.primer_sizes.min, int
        ):
            raise TypeError("Amplicon sizes and primer sizes must be integers")
        if self.gc_clamp[0] > self.gc_clamp[1]:
            raise ValueError("Min primer GC-clamp must be <= max primer GC-clamp")

    def to_input_tags(self) -> dict[Primer3InputTag, Any]:
        """Converts input params to Primer3InputTag to feed directly into Primer3."""
        mapped_dict = {
            Primer3InputTag.PRIMER_NUM_RETURN: self.number_primers_return,
            Primer3InputTag.PRIMER_PRODUCT_OPT_SIZE: self.amplicon_sizes.opt,
            Primer3InputTag.PRIMER_PRODUCT_SIZE_RANGE: (
                f"{self.amplicon_sizes.min}-{self.amplicon_sizes.max}"
            ),
            Primer3InputTag.PRIMER_PRODUCT_MIN_TM: self.amplicon_tms.min,
            Primer3InputTag.PRIMER_PRODUCT_OPT_TM: self.amplicon_tms.opt,
            Primer3InputTag.PRIMER_PRODUCT_MAX_TM: self.amplicon_tms.max,
            Primer3InputTag.PRIMER_MIN_SIZE: self.primer_sizes.min,
            Primer3InputTag.PRIMER_OPT_SIZE: self.primer_sizes.opt,
            Primer3InputTag.PRIMER_MAX_SIZE: self.primer_sizes.max,
            Primer3InputTag.PRIMER_MIN_TM: self.primer_tms.min,
            Primer3InputTag.PRIMER_OPT_TM: self.primer_tms.opt,
            Primer3InputTag.PRIMER_MAX_TM: self.primer_tms.max,
            Primer3InputTag.PRIMER_MIN_GC: self.primer_gcs.min,
            Primer3InputTag.PRIMER_OPT_GC_PERCENT: self.primer_gcs.opt,
            Primer3InputTag.PRIMER_MAX_GC: self.primer_gcs.max,
            Primer3InputTag.PRIMER_GC_CLAMP: self.gc_clamp[0],
            Primer3InputTag.PRIMER_MAX_END_GC: self.gc_clamp[1],
            Primer3InputTag.PRIMER_MAX_POLY_X: self.primer_max_polyX,
            Primer3InputTag.PRIMER_MAX_NS_ACCEPTED: self.primer_max_Ns,
            Primer3InputTag.PRIMER_LOWERCASE_MASKING: 1 if self.avoid_masked_bases else 0,
        }
        return mapped_dict

    @property
    def max_amplicon_length(self) -> int:
        """Max amplicon length"""
        return int(self.amplicon_sizes.max)

    @property
    def max_primer_length(self) -> int:
        """Max primer length"""
        return int(self.primer_sizes.max)

    @property
    def min_primer_length(self) -> int:
        """Minimum primer length."""
        return int(self.primer_sizes.min)
