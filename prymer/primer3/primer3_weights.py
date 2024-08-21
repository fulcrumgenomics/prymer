"""
# Primer3Weights Class and Methods

The Primer3Weights class holds the penalty weights that Primer3 uses to score primer designs.

Primer3 considers the differential between user input (e.g., constraining the optimal
primer size to be 18 bp) and the characteristics of a specific primer design (e.g., if the primer
size is 19 bp). Depending on the "weight" of that characteristic, Primer3 uses an objective function
to score a primer design and help define what an "optimal" design looks like.

By modifying these weights, users can prioritize specific primer design characteristics. Each of
the defaults provided here are derived from the Primer3 manual: https://primer3.org/manual.html

## Examples of interacting with the `Primer3Weights` class


```python
>>> Primer3Weights(product_size_lt=1, product_size_gt=1)
Primer3Weights(product_size_lt=1, product_size_gt=1, ...)
>>> Primer3Weights(product_size_lt=5, product_size_gt=1)
Primer3Weights(product_size_lt=5, product_size_gt=1, ...)

```
"""

from dataclasses import dataclass
from typing import Any

from prymer.primer3.primer3_input_tag import Primer3InputTag


@dataclass(frozen=True, init=True, slots=True)
class Primer3Weights:
    """Holds the weights that Primer3 uses to adjust penalties
     that originate from the designed primer(s).

    The weights that Primer3 uses when a parameter is less than optimal are labeled with "_lt".
    "_gt" weights are penalties applied when a parameter is greater than optimal.

    Please see the Primer3 manual for additional details:
     https://primer3.org/manual.html#globalTags

     Example:
         >>> Primer3Weights() #default implementation
         Primer3Weights(product_size_lt=1, product_size_gt=1, product_tm_lt=0.0, product_tm_gt=0.0, primer_end_stability=0.25, primer_gc_lt=0.25, primer_gc_gt=0.25, primer_self_any=0.1, primer_self_end=0.1, primer_size_lt=0.5, primer_size_gt=0.1, primer_tm_lt=1.0, primer_tm_gt=1.0)

         >>> Primer3Weights(product_size_lt=5)
         Primer3Weights(product_size_lt=5, product_size_gt=1, product_tm_lt=0.0, product_tm_gt=0.0, primer_end_stability=0.25, primer_gc_lt=0.25, primer_gc_gt=0.25, primer_self_any=0.1, primer_self_end=0.1, primer_size_lt=0.5, primer_size_gt=0.1, primer_tm_lt=1.0, primer_tm_gt=1.0)
    """  # noqa: E501

    product_size_lt: int = 1
    product_size_gt: int = 1
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
