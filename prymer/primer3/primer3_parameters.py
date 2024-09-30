"""
# PrimerAndAmpliconParameters and ProbeParameters: Classes and Methods

The [`PrimerAndAmpliconParameters`][prymer.primer3.primer3_parameters.PrimerAndAmpliconParameters]
class stores user input for primer design and maps it to the correct Primer3 fields.

Primer3 considers many criteria for primer design, including characteristics of candidate primers
and the resultant amplicon product, as well as potential complications (off-target priming,
primer dimer formation). Users can specify many of these constraints in Primer3,
some of which are used to quantify a "score" for each primer design.

The PrimerAndAmpliconParameters class stores commonly used constraints for primer design:
GC content, melting temperature, and size of both primers and expected amplicon.
Additional criteria include the maximum homopolymer length, ambiguous bases, and bases in a
dinucleotide run within a primer. By default, primer design avoids masked bases, returns 5 primers,
and sets the GC clamp to be no larger than 5. The PrimerAndAmpliconParameters also stores common
default settings to minimize the tendendancy of primers and probes to anneal to one another. This
self-complementarity can make PCR reactions less efficient and potentially yield nonspecific
amplification. Primer3 supports thermodynamic-based thresholds as well as parameters for maximal
alignment scores.

The `to_input_tags()` method in `PrimerAndAmpliconParameters` converts these parameters into
tag-values pairs for use when executing `Primer3`.

The [`ProbeParameters`][prymer.primer3.primer3_parameters.ProbeParameters]
class stores user input for internal probe design and maps it to the correct Primer3 fields.

Similar to the PrimerAndAmpliconParameters class, the ProbeParameters class can be used to
specify the acceptable ranges of probe sizes, melting temperatures, and GC content. A region can be
excluded from internal probe design based on its start and the length of the region to exclude. This
attribute can help avoid regions that are problematic for oligo design, like low-complexity
sequence tracts.

## Examples

```python
>>> params = PrimerAndAmpliconParameters( \
    amplicon_sizes=MinOptMax(min=100, max=250, opt=200), \
    amplicon_tms=MinOptMax(min=55.0, max=100.0, opt=70.0), \
    primer_sizes=MinOptMax(min=29, max=31, opt=30), \
    primer_tms=MinOptMax(min=63.0, max=67.0, opt=65.0), \
    primer_gcs=MinOptMax(min=30.0, max=65.0, opt=45.0), \
)
>>> for tag, value in params.to_input_tags().items(): \
    print(f"{tag.value} -> {value}")
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
PRIMER_NUM_RETURN -> 5
PRIMER_MAX_SELF_ANY -> 8.0
PRIMER_MAX_SELF_ANY_TH -> 53.0
PRIMER_MAX_SELF_END -> 3.0
PRIMER_MAX_SELF_END_TH -> 53.0
PRIMER_MAX_HAIRPIN_TH -> 53.0

```
"""

import warnings
from dataclasses import dataclass
from dataclasses import field
from typing import Any

from prymer.api.minoptmax import MinOptMax
from prymer.primer3.primer3_input_tag import Primer3InputTag


@dataclass(frozen=True, init=True, slots=True)
class PrimerAndAmpliconParameters:
    """Holds common primer and amplicon design options that Primer3 uses to inform primer design.

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
        primer_max_self_any: the maximal local alignment score of aligning the primer to itself
        primer_max_self_any_thermo: the maximal melting temperature of the most stable structure
            resulting from aligning the primer to itself
        primer_max_self_end: the maximal 3' anchored global alignment score of aligning the primer
            to itself
        primer_max_self_end_thermo: the maximal melting temperature of the most stable structure
            resulting from aligning the 3' end of the primer
        primer_max_hairpin_thermo: the maximal melting temperature of the most stable hairpin
            structure of the primer

    Primer3 uses both thermodynamic and alignment-based approaches to quantify primer
    self-complementarity.

    `primer_max_self_any`, `primer_max_self_any_thermo`, `primer_max_self_end`,
    `primer_max_self_end_thermo`, and `primer_max_hairpin_thermo` are all set to default values as
    specified in the Primer3 manual. The default values of the thermodynamic attributes
    (ending in `_th`) are set to 10 degrees less than the minimal melting temperature specified for
    primer design.

    For `primer_max_self_any` and `primer_max_self_end`, a score of 0.00 indicates that there is no
    reasonable local alignment across the individual primer under consideration.

    In general, these settings are meant to limit problematic oligo self-complementarity
    and avoid primer-dimers or other nonspecific binding of oligos to target sequences.

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
    primer_max_self_any: float = 8.00
    primer_max_self_any_thermo: float = field(init=False)
    primer_max_self_end: float = 3.0
    primer_max_self_end_thermo: float = field(init=False)
    primer_max_hairpin_thermo: float = field(init=False)

    def __post_init__(self) -> None:
        if self.primer_max_dinuc_bases % 2 == 1:
            raise ValueError("Primer Max Dinuc Bases must be an even number of bases")
        if not isinstance(self.amplicon_sizes.min, int) or not isinstance(
            self.primer_sizes.min, int
        ):
            raise TypeError("Amplicon sizes and primer sizes must be integers")
        if self.gc_clamp[0] > self.gc_clamp[1]:
            raise ValueError("Min primer GC-clamp must be <= max primer GC-clamp")

        # Set melting temperature thresholds to be 10 degrees less than the minimum primer tm
        default_thermo_tm: float = self.primer_tms.min - 10.0
        object.__setattr__(self, "primer_max_self_any_thermo", default_thermo_tm)
        object.__setattr__(self, "primer_max_self_end_thermo", default_thermo_tm)
        object.__setattr__(self, "primer_max_hairpin_thermo", default_thermo_tm)

    def to_input_tags(self) -> dict[Primer3InputTag, Any]:
        """Converts input params to Primer3InputTag to feed directly into Primer3."""
        mapped_dict: dict[Primer3InputTag, Any] = {
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
            Primer3InputTag.PRIMER_NUM_RETURN: self.number_primers_return,
            Primer3InputTag.PRIMER_MAX_SELF_ANY: self.primer_max_self_any,
            Primer3InputTag.PRIMER_MAX_SELF_ANY_TH: self.primer_max_self_any_thermo,
            Primer3InputTag.PRIMER_MAX_SELF_END: self.primer_max_self_end,
            Primer3InputTag.PRIMER_MAX_SELF_END_TH: self.primer_max_self_end_thermo,
            Primer3InputTag.PRIMER_MAX_HAIRPIN_TH: self.primer_max_hairpin_thermo,
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


@dataclass(frozen=True, init=True, slots=True)
class Primer3Parameters(PrimerAndAmpliconParameters):
    """A deprecated alias for `PrimerAndAmpliconParameters` intended to maintain backwards
    compatibility with earlier releases of `prymer`."""

    warnings.warn(
        "The Primer3Parameters class was deprecated, use PrimerAndAmpliconParameters instead",
        DeprecationWarning,
        stacklevel=2,
    )


@dataclass(frozen=True, init=True, slots=True)
class ProbeParameters:
    """Holds common primer design options that Primer3 uses to inform internal probe design.

    Attributes:
        probe_sizes: the min, optimal, and max probe size
        probe_tms: the min, optimal, and max probe melting temperatures
        probe_gcs: the min and max GC content for individual probes
        number_probes_return: the number of probes to return
        probe_max_dinuc_bases: the max  number of bases in a dinucleotide run in a probe
        probe_max_polyX: the max homopolymer length acceptable within a probe
        probe_max_Ns: the max number of ambiguous bases acceptable within a probe
        probe_max_self_any: the maximal local alignment score of aligning the probe to itself
        probe_max_self_any_thermo: the maximal melting temperature of the most stable structure
            resulting from aligning the probe to itself
        probe_max_self_end: max allowable 3'-anchored global alignment score when testing a
            probe for self-complementarity
        probe_max_self_end_thermo: the maximal melting temperature of the most stable structure
            resulting from aligning the 3' end of the probe
        probe_max_hairpin_thermo: most stable monomer structure as calculated by a thermodynamic
            approach

    The attributes that have default values specified take their default values from the
    Primer3 manual.

    Please see the Primer3 manual for additional details: https://primer3.org/manual.html#globalTags

    Primer3 uses both thermodynamic and alignment-based approaches to quantify oligo
    self-complementarity.

    `primer_max_self_any`, `primer_max_self_any_thermo`, `primer_max_self_end`,
    `primer_max_self_end_thermo`, and `primer_max_hairpin_thermo` are all set to default values as
    specified in the Primer3 manual. The default values of the thermodynamic attributes
    (ending in `_th`) are set to 10 degrees less than the minimal melting temperature specified for
    primer design.

    For `probe_max_self_any` and `probe_max_self_end`, a score of 0.00 indicates that there is no
    reasonable local alignment across the individual primer under consideration. These scores are
    always positive.

    In general, these settings are meant to limit problematic oligo self-complementarity
    and avoid primer-dimers or other nonspecific binding of probes to target sequences.

    Note that the Primer3 documentation advises that, while `probe_max_end_any` is meaningless
    when applied to internal probes used for hybridization-based detection,
    `PRIMER_INTERNAL_MAX_SELF_END` should be set at least as high as `PRIMER_INTERNAL_MAX_SELF_ANY`.
    Therefore, both parameters are exposed here.

    """

    probe_sizes: MinOptMax[int]
    probe_tms: MinOptMax[float]
    probe_gcs: MinOptMax[float]
    number_probes_return: int = 5
    probe_max_dinuc_bases: int = 4
    probe_max_polyX: int = 5
    probe_max_Ns: int = 0
    probe_max_self_any: float = 12.0
    probe_max_self_any_thermo: float = field(init=False)
    probe_max_self_end: float = 12.0
    probe_max_self_end_thermo: float = field(init=False)
    probe_max_hairpin_thermo: float = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.probe_sizes.min, int):
            raise TypeError("Probe sizes must be integers")
        if not isinstance(self.probe_tms.min, float) or not isinstance(self.probe_gcs.min, float):
            raise TypeError("Probe melting temperatures and GC content must be floats")
        if self.probe_max_dinuc_bases % 2 == 1:
            raise ValueError("Max threshold for dinucleotide bases must be an even number of bases")

        # Set melting temperature thresholds to be 10 degrees less than the minimum primer tm
        default_thermo_tm: float = self.probe_tms.min - 10.0
        object.__setattr__(self, "probe_max_self_any_thermo", default_thermo_tm)
        object.__setattr__(self, "probe_max_self_end_thermo", default_thermo_tm)
        object.__setattr__(self, "probe_max_hairpin_thermo", default_thermo_tm)

    def to_input_tags(self) -> dict[Primer3InputTag, Any]:
        """Converts input params to Primer3InputTag to feed directly into Primer3."""
        mapped_dict: dict[Primer3InputTag, Any] = {
            Primer3InputTag.PRIMER_INTERNAL_MIN_SIZE: self.probe_sizes.min,
            Primer3InputTag.PRIMER_INTERNAL_OPT_SIZE: self.probe_sizes.opt,
            Primer3InputTag.PRIMER_INTERNAL_MAX_SIZE: self.probe_sizes.max,
            Primer3InputTag.PRIMER_INTERNAL_MIN_TM: self.probe_tms.min,
            Primer3InputTag.PRIMER_INTERNAL_OPT_TM: self.probe_tms.opt,
            Primer3InputTag.PRIMER_INTERNAL_MAX_TM: self.probe_tms.max,
            Primer3InputTag.PRIMER_INTERNAL_MIN_GC: self.probe_gcs.min,
            Primer3InputTag.PRIMER_INTERNAL_OPT_GC_PERCENT: self.probe_gcs.opt,
            Primer3InputTag.PRIMER_INTERNAL_MAX_GC: self.probe_gcs.max,
            Primer3InputTag.PRIMER_INTERNAL_MAX_POLY_X: self.probe_max_polyX,
            Primer3InputTag.PRIMER_INTERNAL_MAX_NS_ACCEPTED: self.probe_max_Ns,
            Primer3InputTag.PRIMER_INTERNAL_MAX_SELF_ANY: self.probe_max_self_any,
            Primer3InputTag.PRIMER_INTERNAL_MAX_SELF_ANY_TH: self.probe_max_self_any_thermo,
            Primer3InputTag.PRIMER_INTERNAL_MAX_SELF_END: self.probe_max_self_end,
            Primer3InputTag.PRIMER_INTERNAL_MAX_SELF_END_TH: self.probe_max_self_end_thermo,
            Primer3InputTag.PRIMER_INTERNAL_MAX_HAIRPIN_TH: self.probe_max_hairpin_thermo,
        }

        return mapped_dict
