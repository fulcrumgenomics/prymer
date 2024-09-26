"""
# Primer3Input Class and Methods

This module contains the [`Primer3Input`][prymer.primer3.Primer3Input] class. The class
wraps together different helper classes to assemble user-specified criteria and parameters for
input to Primer3.

The module uses:

1. [`PrimerAndAmpliconParameters`][prymer.primer3.primer3_parameters.Primer3Parameters]
    to specify user-specified criteria for primer design
2. [`ProbeParameters`][prymer.primer3.primer3_parameters.ProbeParameters]
    to specify user-specified criteria for probe design
3. [`PrimerAndAmpliconWeights`][prymer.primer3.primer3_weights.PrimerAndAmpliconWeights]
    to establish penalties based on those criteria
4. [`ProbeWeights`][prymer.primer3.primer3_weights.ProbeWeights] to specify penalties based on probe
    design criteria
5. [`Primer3Task`][prymer.primer3.primer3_task.Primer3Task] to organize task-specific
    logic.
6. [`Span`](index.md#prymer.api.span.Span] to specify the target region.

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
>>> params = PrimerAndAmpliconParameters( \
    amplicon_sizes=MinOptMax(min=100, max=250, opt=200), \
    amplicon_tms=MinOptMax(min=55.0, max=100.0, opt=70.0), \
    primer_sizes=MinOptMax(min=29, max=31, opt=30), \
    primer_tms=MinOptMax(min=63.0, max=67.0, opt=65.0), \
    primer_gcs=MinOptMax(min=30.0, max=65.0, opt=45.0), \
    )
>>> design_input = Primer3Input(target=target, \
                            primer_and_amplicon_params=params, \
                            task=DesignLeftPrimersTask() \
    )

>>> for tag, value in design_input.to_input_tags(design_region=design_region).items(): \
    print(f"{tag.value} -> {value}")
PRIMER_TASK -> pick_primer_list
PRIMER_PICK_LEFT_PRIMER -> 1
PRIMER_PICK_RIGHT_PRIMER -> 0
PRIMER_PICK_INTERNAL_OLIGO -> 0
SEQUENCE_INCLUDED_REGION -> 1,51
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
PRIMER_WT_SELF_ANY_TH -> 0.0
PRIMER_WT_SELF_END_TH -> 0.0
PRIMER_WT_HAIRPIN_TH -> 0.0
"""

from dataclasses import MISSING
from dataclasses import dataclass
from dataclasses import fields
from typing import Any
from typing import Optional

from prymer.api.span import Span
from prymer.primer3.primer3_input_tag import Primer3InputTag
from prymer.primer3.primer3_parameters import PrimerAndAmpliconParameters
from prymer.primer3.primer3_parameters import ProbeParameters
from prymer.primer3.primer3_task import Primer3TaskType
from prymer.primer3.primer3_weights import PrimerAndAmpliconWeights
from prymer.primer3.primer3_weights import ProbeWeights


@dataclass(frozen=True, init=True, slots=True)
class Primer3Input:
    """Assembles necessary inputs for Primer3 to orchestrate primer, primer pair, and/or internal
    probe design."""

    target: Span
    task: Primer3TaskType
    primer_and_amplicon_params: Optional[PrimerAndAmpliconParameters] = None
    probe_params: Optional[ProbeParameters] = None
    primer_weights: Optional[PrimerAndAmpliconWeights] = None
    probe_weights: Optional[ProbeWeights] = None

    def __post_init__(self) -> None:
        # check for at least one set of params
        # for the set of params given, check that weights were given; use defaults if not given
        if self.primer_and_amplicon_params is None and self.probe_params is None:
            raise ValueError(
                "Primer3 requires at least one set of parameters"
                " for either primer or probe design"
            )

        if self.primer_and_amplicon_params is not None and self.primer_weights is None:
            object.__setattr__(self, "primer_weights", PrimerAndAmpliconWeights())

        if self.probe_params is not None and self.probe_weights is None:
            object.__setattr__(self, "probe_weights", ProbeWeights())

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
        assembled_tags: dict[Primer3InputTag, Any] = {**primer3_task_params}

        optional_attributes = {
            field.name: getattr(self, field.name)
            for field in fields(self)
            if field.default is not MISSING
        }
        for settings in optional_attributes.values():
            if settings is not None:
                assembled_tags.update(settings.to_input_tags())

        return assembled_tags
