"""
# Primer3Input Class and Methods

This module contains the [`Primer3Input`][prymer.primer3.Primer3Input] class. The class
wraps together different helper classes to assemble user-specified criteria and parameters for
input to Primer3.

The module uses:

1. [`Primer3Parameters`][prymer.primer3.primer3_parameters.Primer3Parameters]
to specify user-specified criteria for primer design
2. [`Primer3Weights`][prymer.primer3.primer3_weights.Primer3Weights] to establish penalties
based on those criteria
3. [`Primer3Task`][prymer.primer3.primer3_task.Primer3Task] to organize task-specific
    logic.
4. [`Span`](index.md#prymer.api.span.Span] to specify the target region.

The `Primer3Input.to_input_tags(]` method
The main purpose of this class is to generate the
[`Primer3InputTag`s][prymer.primer3.primer3_input_tag.Primer3InputTag]s required by
`Primer3` for specifying how to design the primers, returned by the `to_input_tags(]` method.

## Examples

The following examples builds the `Primer3` tags for designing left primers:

```python
>>> from prymer.api import MinOptMax, Strand
>>> from prymer.primer3 import DesignLeftPrimersTask
>>> target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)
>>> design_region = Span(refname="chr1", start=150, end=300, strand=Strand.POSITIVE)
>>> params = Primer3Parameters( \
    amplicon_sizes=MinOptMax(min=100, max=250, opt=200), \
    amplicon_tms=MinOptMax(min=55.0, max=100.0, opt=70.0), \
    primer_sizes=MinOptMax(min=29, max=31, opt=30), \
    primer_tms=MinOptMax(min=63.0, max=67.0, opt=65.0), \
    primer_gcs=MinOptMax(min=30.0, max=65.0, opt=45.0), \
)
>>> design_input = Primer3Input(target=target, params=params, task=DesignLeftPrimersTask())
>>> for tag, value in design_input.to_input_tags(design_region=design_region).items(): \
    print(f"{tag.value} -> {value}")
PRIMER_TASK -> pick_primer_list
PRIMER_PICK_LEFT_PRIMER -> 1
PRIMER_PICK_RIGHT_PRIMER -> 0
PRIMER_PICK_INTERNAL_OLIGO -> 0
SEQUENCE_INCLUDED_REGION -> 1,51
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
PRIMER_PAIR_WT_PRODUCT_SIZE_LT -> 1
PRIMER_PAIR_WT_PRODUCT_SIZE_GT -> 1
PRIMER_PAIR_WT_PRODUCT_TM_LT -> 0.0
PRIMER_PAIR_WT_PRODUCT_TM_GT -> 0.0
PRIMER_WT_END_STABILITY -> 0.25
PRIMER_WT_GC_PERCENT_LT -> 0.25
PRIMER_WT_GC_PERCENT_GT -> 0.25
PRIMER_WT_SELF_ANY -> 0.1
PRIMER_WT_SELF_END -> 0.1
PRIMER_WT_SIZE_LT -> 0.5
PRIMER_WT_SIZE_GT -> 0.1
PRIMER_WT_TM_LT -> 1.0
PRIMER_WT_TM_GT -> 1.0
"""

from dataclasses import dataclass
from typing import Any

from prymer.api.span import Span
from prymer.primer3.primer3_input_tag import Primer3InputTag
from prymer.primer3.primer3_parameters import Primer3Parameters
from prymer.primer3.primer3_task import Primer3TaskType
from prymer.primer3.primer3_weights import Primer3Weights


@dataclass(frozen=True, init=True, slots=True)
class Primer3Input:
    """Assembles necessary inputs for Primer3 to orchestrate primer and/or primer pair design."""

    target: Span
    task: Primer3TaskType
    params: Primer3Parameters
    weights: Primer3Weights = Primer3Weights()

    def to_input_tags(self, design_region: Span) -> dict[Primer3InputTag, Any]:
        """Assembles `Primer3InputTag` and values for input to `Primer3`

        The target region must be wholly contained within design region.

        Args:
            design_region: the design region, which wholly contains the target region, in which
                    primers are to be designed.

        Returns:
            a mapping of `Primer3InputTag`s to associated value
        """
        primer3_task_params = self.task.to_input_tags(
            design_region=design_region, target=self.target
        )
        assembled_tags = {
            **primer3_task_params,
            **self.params.to_input_tags(),
            **self.weights.to_input_tags(),
        }
        return assembled_tags
