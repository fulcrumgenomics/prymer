"""
# MinOptMax Classes and Methods

This module contains a class and class methods to hold a range of user-specific thresholds.

Each of the three class attributes represents a minimal, optimal, and maximal value.
The three values can be either int or float values but all must be of the same type within one
MinOptMax object (for example, `min` cannot be a float while `max` is an int).

Primer3 will use these values downstream to set an allowable range of specific parameters that
inform primer design. For example, Primer3 can constrain primer melting temperature
to be within a range of 55.0 - 65.0 Celsius based on an input
`MinOptMax(min=55.0, opt=60.0, max=65.0)`.

## Examples of interacting with the `MinOptMax` class

```python
>>> thresholds = MinOptMax(min=1.0, opt=2.0, max=4.0)
>>> print(thresholds)
(min:1.0, opt:2.0, max:4.0)
>>> list(thresholds)
[1.0, 2.0, 4.0]

```
"""

from dataclasses import dataclass
from typing import Generic
from typing import Iterator
from typing import TypeVar

Numeric = TypeVar("Numeric", int, float)


@dataclass(slots=True, frozen=True, init=True)
class MinOptMax(Generic[Numeric]):
    """Stores a minimum, an optimal, and a maximum value (either all ints or all floats).

    `min` must be less than `max`. `opt` should be greater than the `min`
     value and less than the `max` value.

    Attributes:
        min: the minimum value
        opt: the optimal value
        max: the maximum value


    Raises:
        ValueError: if min > max
        ValueError: if `min` is not less than `opt` and `opt` is not less than `max`
    """

    min: Numeric
    opt: Numeric
    max: Numeric

    def __post_init__(self) -> None:
        dtype = type(self.min)
        if not isinstance(self.max, dtype) or not isinstance(self.opt, dtype):
            raise TypeError(
                "Min, opt, and max must be the same type; "
                f"received min: {dtype}, opt: {type(self.opt)}, max: {type(self.max)}"
            )
        if self.min > self.max:
            raise ValueError(
                f"Min must be no greater than max; received min: {self.min}, max: {self.max}"
            )
        if not (self.min <= self.opt <= self.max):
            raise ValueError(
                "Arguments must satisfy min <= opt <= max: "
                f"received min: {self.min}, opt: {self.opt}, max: {self.max}"
            )

    def __iter__(self) -> Iterator[float]:
        """Returns an iterator of min, opt, and max"""
        return iter([self.min, self.opt, self.max])

    def __str__(self) -> str:
        """Returns a string representation of min, opt, and max"""
        return f"(min:{self.min}, opt:{self.opt}, max:{self.max})"
