"""
# Primer3Weights Class and Methods

The PrimerAndAmpliconWeights class holds the penalty weights that Primer3 uses to score
primer designs.

The ProbeWeights class holds the penalty weights that Primer3 uses to score internal probe designs.

Primer3 considers the differential between user input (e.g., constraining the optimal
primer size to be 18 bp) and the characteristics of a specific primer design (e.g., if the primer
size is 19 bp). Depending on the "weight" of that characteristic, Primer3 uses an objective function
to score a primer design and help define what an "optimal" design looks like.

By modifying these weights, users can prioritize specific primer design characteristics. Each of
the defaults provided here are derived from the Primer3 manual: https://primer3.org/manual.html

## Examples of interacting with the `PrimerAndAmpliconWeights` class

Example:
>>> PrimerAndAmpliconWeights() # default implementation
PrimerAndAmpliconWeights(product_size_lt=1.0, product_size_gt=1, product_tm_lt=0.0, product_tm_gt=0.0, primer_end_stability=0.25, primer_gc_lt=0.25, primer_gc_gt=0.25, primer_self_any=0.1, primer_self_end=0.1, primer_size_lt=0.5, primer_size_gt=0.1, primer_tm_lt=1.0, primer_tm_gt=1.0)
>>> PrimerAndAmpliconWeights(product_size_lt=5)
PrimerAndAmpliconWeights(product_size_lt=5.0, product_size_gt=1, product_tm_lt=0.0, product_tm_gt=0.0, primer_end_stability=0.25, primer_gc_lt=0.25, primer_gc_gt=0.25, primer_self_any=0.1, primer_self_end=0.1, primer_size_lt=0.5, primer_size_gt=0.1, primer_tm_lt=1.0, primer_tm_gt=1.0)

"""  # noqa: E501

from dataclasses import dataclass
from typing import Any

from prymer.primer3.primer3_input_tag import Primer3InputTag


@dataclass(frozen=True, init=True, slots=True)
class PrimerAndAmpliconWeights:
    """Holds the primer-specific weights that Primer3 uses to adjust design penalties.

    The weights that Primer3 uses when a parameter is less than optimal are labeled with "_lt".
    "_gt" weights are penalties applied when a parameter is greater than optimal.

    Some of these settings depart from the default settings enumerated in the Primer3 manual.
    Please see the Primer3 manual for additional details:
     https://primer3.org/manual.html#globalTags

    Attributes:
        product_size_lt: weight for products shorter than
            `PrimerAndAmpliconParameters.amplicon_sizes.opt`
        product_size_gt: weight for products longer than
            `PrimerAndAmpliconParameters.amplicon_sizes.opt`
        product_tm_lt: weight for products with a Tm lower than
            `PrimerAndAmpliconParameters.amplicon_tms.opt`
        product_tm_gt: weight for products with a Tm greater than
            `PrimerAndAmpliconParameters.amplicon_tms.opt`
        primer_end_stability: penalty for the calculated maximum stability
            for the last five 3' bases of primer
        primer_gc_lt: penalty for primers with GC percent lower than
            `PrimerAndAmpliconParameters.primer_gcs.opt`
        primer_gc_gt: penalty weight for primers with GC percent higher than
            `PrimerAndAmpliconParameters.primer_gcs.opt`
     Example:
         >>> PrimerAndAmpliconWeights() #default implementation
         Primer3Weights(product_size_lt=1.0, product_size_gt=1.0, product_tm_lt=0.0, product_tm_gt=0.0, primer_end_stability=0.25, primer_gc_lt=0.25, primer_gc_gt=0.25, primer_self_any=0.1, primer_self_end=0.1, primer_size_lt=0.5, primer_size_gt=0.1, primer_tm_lt=1.0, primer_tm_gt=1.0)
         >>> PrimerAndAmpliconWeights(product_size_lt=5.0)
         Primer3Weights(product_size_lt=5.0, product_size_gt=1.0, product_tm_lt=0.0, product_tm_gt=0.0, primer_end_stability=0.25, primer_gc_lt=0.25, primer_gc_gt=0.25, primer_self_any=0.1, primer_self_end=0.1, primer_size_lt=0.5, primer_size_gt=0.1, primer_tm_lt=1.0, primer_tm_gt=1.0)
    """  # noqa: E501

    product_size_lt: float = 1.0
    product_size_gt: float = 1.0
    product_tm_lt: float = 0.0
    product_tm_gt: float = 0.0
    primer_end_stability: float = 0.25
    primer_gc_lt: float = 0.25
    primer_gc_gt: float = 0.25
    primer_self_any: float = 0.1
    primer_self_end: float = 0.1
    primer_size_lt: float = 0.5
    primer_size_gt: float = 0.1
    primer_tm_lt: float = 1.0
    primer_tm_gt: float = 1.0

    def to_input_tags(self) -> dict[Primer3InputTag, Any]:
        """Maps weights to Primer3InputTag to feed directly into Primer3."""
        mapped_dict = {
            Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_LT: self.product_size_lt,
            Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_GT: self.product_size_gt,
            Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_TM_LT: self.product_tm_lt,
            Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_TM_GT: self.product_tm_gt,
            Primer3InputTag.PRIMER_WT_END_STABILITY: self.primer_end_stability,
            Primer3InputTag.PRIMER_WT_GC_PERCENT_LT: self.primer_gc_lt,
            Primer3InputTag.PRIMER_WT_GC_PERCENT_GT: self.primer_gc_gt,
            Primer3InputTag.PRIMER_WT_SELF_ANY: self.primer_self_any,
            Primer3InputTag.PRIMER_WT_SELF_END: self.primer_self_end,
            Primer3InputTag.PRIMER_WT_SIZE_LT: self.primer_size_lt,
            Primer3InputTag.PRIMER_WT_SIZE_GT: self.primer_size_gt,
            Primer3InputTag.PRIMER_WT_TM_LT: self.primer_tm_lt,
            Primer3InputTag.PRIMER_WT_TM_GT: self.primer_tm_gt,
        }
        return mapped_dict


@dataclass(frozen=True, init=True, slots=True)
class ProbeWeights:
    """Holds the probe-specific weights that Primer3 uses to adjust design penalties.

    Attributes:
        probe_size_lt: penalty for probes shorter than `ProbeParameters.probe_sizes.opt`
        probe_size_gt: penalty for probes longer than `ProbeParameters.probe_sizes.opt`
        probe_tm_lt: penalty for probes with a Tm lower than `ProbeParameters.probe_tms.opt`
        probe_tm_gt: penalty for probes with a Tm greater than `ProbeParameters.probe_tms.opt`
        probe_gc_lt: penalty for probes with GC content lower than `ProbeParameters.probe_gcs.opt`
        probe_gc_gt: penalty for probes with GC content greater than `ProbeParameters.probe_gcs.opt`

    """

    probe_size_lt: float = 0.25
    probe_size_gt: float = 0.25
    probe_tm_lt: float = 1.0
    probe_tm_gt: float = 1.0
    probe_gc_lt: float = 0.5
    probe_gc_gt: float = 0.5

    def to_input_tags(self) -> dict[Primer3InputTag, Any]:
        """Maps weights to Primer3InputTag to feed directly into Primer3."""
        mapped_dict = {
            Primer3InputTag.PRIMER_INTERNAL_WT_SIZE_LT: self.probe_size_lt,
            Primer3InputTag.PRIMER_INTERNAL_WT_SIZE_GT: self.probe_size_gt,
            Primer3InputTag.PRIMER_INTERNAL_WT_TM_LT: self.probe_tm_lt,
            Primer3InputTag.PRIMER_INTERNAL_WT_TM_GT: self.probe_tm_gt,
            Primer3InputTag.PRIMER_INTERNAL_WT_GC_PERCENT_LT: self.probe_gc_lt,
            Primer3InputTag.PRIMER_INTERNAL_WT_GC_PERCENT_GT: self.probe_gc_gt,
        }
        return mapped_dict
