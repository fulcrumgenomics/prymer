"""
# Primer3Task Class and Methods

This module contains a Primer3Task class to provide design-specific parameters to Primer3.  The
classes are primarily used in [`Primer3Input`][prymer.primer3.Primer3Input].

Primer3 can design single primers ("left" and "right") as well as primer pairs.
The design task "type" dictates which type of primers to pick and informs the design region.
These parameters are aligned to the correct Primer3 settings and fed directly into Primer3.

Three types of tasks are available:

1. [`DesignPrimerPairsTask`][prymer.primer3.primer3_task.DesignPrimerPairsTask] -- task
    for designing _primer pairs_.
2. [`DesignLeftPrimersTask`][prymer.primer3.primer3_task.DesignLeftPrimersTask] -- task
    for designing primers to the _left_ (5') of the design region on the top/positive strand.
3. [`DesignRightPrimersTask`][prymer.primer3.primer3_task.DesignRightPrimersTask] -- task
    for designing primers to the _right_ (3') of the design region on the bottom/negative strand.

The main purpose of these classes are to generate the
[`Primer3InputTag`s][prymer.primer3.primer3_input_tag.Primer3InputTag]s required by
`Primer3` for specifying how to design the primers, returned by the `to_input_tags()` method.  The
target and design region are provided to this method, where the target region is wholly contained
within design region.  This leaves a left and right primer region respectively, that are
the two regions that remain after removing the wholly contained (inner) target regions.

Therefore, the left primers shall be designed from the start of the design region to the
start of the target region, while right primers shall be designed from the end of the target region
through to the end of the design region.

In addition to the `to_input_tags()` method, each `Primer3TaskType` provides a `task_type` and
`count_tag` class property.  The former is a `TaskType` enumeration that represents the type of
design task, while the latter is the tag returned by Primer3 that provides the number of primers
returned.

## Examples

Suppose we have the following design and target regions:

```python
>>> from prymer.api import Strand
>>> design_region = Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE)
>>> target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
>>> design_region.contains(target)
True

```

The task for designing primer pairs:

```python
>>> task = DesignPrimerPairsTask()
>>> task.task_type.value
'PAIR'
>>> task.count_tag
'PRIMER_PAIR_NUM_RETURNED'
>>> [(tag.value, value) for tag, value in task.to_input_tags(target=target, design_region=design_region).items()]
[('PRIMER_TASK', 'generic'), ('PRIMER_PICK_LEFT_PRIMER', 1), ('PRIMER_PICK_RIGHT_PRIMER', 1), ('PRIMER_PICK_INTERNAL_OLIGO', 0), ('SEQUENCE_TARGET', '200,101')]

```

The tasks for designing left primers:

```python
>>> task = DesignLeftPrimersTask()
>>> task.task_type.value
'LEFT'
>>> task.count_tag
'PRIMER_LEFT_NUM_RETURNED'
>>> [(tag.value, value) for tag, value in task.to_input_tags(target=target, design_region=design_region).items()]
[('PRIMER_TASK', 'pick_primer_list'), ('PRIMER_PICK_LEFT_PRIMER', 1), ('PRIMER_PICK_RIGHT_PRIMER', 0), ('PRIMER_PICK_INTERNAL_OLIGO', 0), ('SEQUENCE_INCLUDED_REGION', '1,199')]

```

The tasks for designing left primers:

```python
>>> task = DesignRightPrimersTask()
>>> task.task_type.value
'RIGHT'
>>> task.count_tag
'PRIMER_RIGHT_NUM_RETURNED'
>>> [(tag.value, value) for tag, value in task.to_input_tags(target=target, design_region=design_region).items()]
[('PRIMER_TASK', 'pick_primer_list'), ('PRIMER_PICK_LEFT_PRIMER', 0), ('PRIMER_PICK_RIGHT_PRIMER', 1), ('PRIMER_PICK_INTERNAL_OLIGO', 0), ('SEQUENCE_INCLUDED_REGION', '300,200')]

```

"""  # noqa: E501

from abc import ABC
from abc import abstractmethod
from enum import auto
from enum import unique
from typing import Any
from typing import ClassVar
from typing import TypeAlias
from typing import Union
from typing import final

from strenum import UppercaseStrEnum

from prymer.api.span import Span
from prymer.primer3.primer3_input_tag import Primer3InputTag

Primer3TaskType: TypeAlias = Union[
    "DesignPrimerPairsTask", "DesignLeftPrimersTask", "DesignRightPrimersTask"
]
"""Type alias for all `Primer3Task`s, to enable exhaustiveness checking."""


@unique
class TaskType(UppercaseStrEnum):
    """Represents the type of design task, either design primer pairs, or individual primers
    (left or right)."""

    # Developer Note: the names of this enum are important, as they are used as-is for the
    # count_tag in `Primer3Task`.

    PAIR = auto()
    LEFT = auto()
    RIGHT = auto()


class Primer3Task(ABC):
    """Abstract base class from which the other classes derive."""

    @final
    def to_input_tags(self, target: Span, design_region: Span) -> dict[Primer3InputTag, Any]:
        """Returns the input tags specific to this type of task for Primer3.

        Ensures target region is wholly contained within design region. Subclass specific
        implementation aligns the set of input parameters specific to primer pair or
        single primer design.

        This implementation mimics the Template method, a Behavioral Design Pattern in which the
        abstract base class contains a rough skeleton of methods and the derived subclasses
        implement the details of those methods. In this case, each of the derived subclasses
        will first use the base class `to_input_tags()` method to check that the target region is
        wholly contained within the design region. If so, they will implement task-specific logic
        for the `to_input_tags()` method.

        Args:
            target: the target region (to be amplified)
            design_region: the design region, which wholly contains the target region, in which
                primers are to be designed


        Raises:
            ValueError: if the target region is not contained within the design region

        Returns:
            The input tags for Primer3
        """
        if not design_region.contains(target):
            raise ValueError(
                "Target not contained within design region: "
                f"target:{target.__str__()},"
                f"design_region: {design_region.__str__()}"
            )

        return self._to_input_tags(target=target, design_region=design_region)

    task_type: ClassVar[TaskType] = NotImplemented
    """Tracks task type for downstream analysis"""

    count_tag: ClassVar[str] = NotImplemented
    """The tag returned by Primer3 that provides the number of primers returned"""

    @classmethod
    @abstractmethod
    def _to_input_tags(cls, target: Span, design_region: Span) -> dict[Primer3InputTag, Any]:
        """Aligns the set of input parameters specific to primer pair or single primer design"""

    @classmethod
    def __init_subclass__(cls, task_type: TaskType, **kwargs: Any) -> None:
        # See: https://docs.python.org/3/reference/datamodel.html#object.__init_subclass__
        super().__init_subclass__(**kwargs)

        cls.task_type = task_type
        cls.count_tag = f"PRIMER_{task_type}_NUM_RETURNED"


class DesignPrimerPairsTask(Primer3Task, task_type=TaskType.PAIR):
    """Stores task-specific Primer3 settings for designing primer pairs"""

    @classmethod
    def _to_input_tags(cls, target: Span, design_region: Span) -> dict[Primer3InputTag, Any]:
        return {
            Primer3InputTag.PRIMER_TASK: "generic",
            Primer3InputTag.PRIMER_PICK_LEFT_PRIMER: 1,
            Primer3InputTag.PRIMER_PICK_RIGHT_PRIMER: 1,
            Primer3InputTag.PRIMER_PICK_INTERNAL_OLIGO: 0,
            Primer3InputTag.SEQUENCE_TARGET: f"{target.start - design_region.start + 1},"
            f"{target.length}",
        }


class DesignLeftPrimersTask(Primer3Task, task_type=TaskType.LEFT):
    """Stores task-specific characteristics for designing left primers."""

    @classmethod
    def _to_input_tags(cls, target: Span, design_region: Span) -> dict[Primer3InputTag, Any]:
        return {
            Primer3InputTag.PRIMER_TASK: "pick_primer_list",
            Primer3InputTag.PRIMER_PICK_LEFT_PRIMER: 1,
            Primer3InputTag.PRIMER_PICK_RIGHT_PRIMER: 0,
            Primer3InputTag.PRIMER_PICK_INTERNAL_OLIGO: 0,
            Primer3InputTag.SEQUENCE_INCLUDED_REGION: f"1,{target.start - design_region.start}",
        }


class DesignRightPrimersTask(Primer3Task, task_type=TaskType.RIGHT):
    """Stores task-specific characteristics for designing right primers"""

    @classmethod
    def _to_input_tags(cls, target: Span, design_region: Span) -> dict[Primer3InputTag, Any]:
        start = target.end - design_region.start + 1
        length = design_region.end - target.end
        return {
            Primer3InputTag.PRIMER_TASK: "pick_primer_list",
            Primer3InputTag.PRIMER_PICK_LEFT_PRIMER: 0,
            Primer3InputTag.PRIMER_PICK_RIGHT_PRIMER: 1,
            Primer3InputTag.PRIMER_PICK_INTERNAL_OLIGO: 0,
            Primer3InputTag.SEQUENCE_INCLUDED_REGION: f"{start},{length}",
        }
