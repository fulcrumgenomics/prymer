"""
# Class and Methods for primer-like objects

The `PrimerLike` class is an abstract base class designed to represent primer-like objects,
such as individual primers or primer pairs. This class encapsulates common attributes and
provides a foundation for more specialized implementations.

In particular, the following methods/attributes need to be implemented:

- [`span()`][prymer.api.primer_like.PrimerLike.span] -- the mapping of the primer-like
    object to the genome.
- [`bases()`][prymer.api.primer_like.PrimerLike.bases] -- the bases of the primer-like
    object, or `None` if not available.
- [`to_bed12_row()`][prymer.api.primer_like.PrimerLike.to_bed12_row] -- the 12-field BED
    representation of this primer-like object.

See the following concrete implementations:

- [`Primer`][prymer.api.primer.Primer] -- a class to store an individual primer
- [`PrimerPair`][prymer.api.primer_pair.PrimerPair] -- a class to store a primer pair

"""

from abc import ABC
from abc import abstractmethod
from dataclasses import dataclass
from typing import Optional
from typing import TypeVar
from typing import assert_never

from fgpyo.sequence import gc_content

from prymer.api.span import Span
from prymer.api.span import Strand

MISSING_BASES_STRING: str = "*"
"""Used in string representations of primer-like objects when their `bases` property return None"""


@dataclass(frozen=True, init=True, slots=True)
class PrimerLike(ABC):
    """
    An abstract base class for primer-like objects, such as individual primers or primer pairs.

    Attributes:
        name: an optional name to use for the primer

    The `id` field shall be the 'name' field if supplied, or a generated value based on the
    location of the primer-like object.
    """

    name: Optional[str] = None

    def __post_init__(self) -> None:
        # If supplied, bases must be a non-empty String & the same length as the
        # span length
        if self.bases is not None:
            if len(self.bases) == 0:
                raise ValueError("Bases must not be an empty string")
            elif self.span.length != len(self.bases):
                raise ValueError(
                    "Conflicting lengths: "
                    f"span length={self.span.length},"
                    f" sequence length={len(self.bases)}"
                )

    @property
    @abstractmethod
    def span(self) -> Span:
        """Returns the mapping of the primer-like object to a genome."""

    @property
    @abstractmethod
    def bases(self) -> Optional[str]:
        """Returns the base sequence of the primer-like object."""

    @property
    def percent_gc_content(self) -> float:
        """
        The GC of the amplicon sequence in the range 0-100, or zero if there is no amplicon
        sequence.
        """
        if self.bases is None:
            return 0.0
        else:
            return round(gc_content(self.bases) * 100, 3)

    @property
    def id(self) -> str:
        """
        Returns the identifier for the primer-like object. This shall be the `name`
        if one exists, otherwise a generated value based on the location of the object.
        """
        if self.name is not None:
            return self.name
        else:
            return self.location_string

    @property
    def location_string(self) -> str:
        """Returns a string representation of the location of the primer-like object."""
        return (
            f"{self.span.refname}_{self.span.start}_"
            + f"{self.span.end}_{self._strand_to_location_string()}"
        )

    @abstractmethod
    def to_bed12_row(self) -> str:
        """
        Formats the primer-like into 12 tab-separated fields matching the BED 12-column spec.
        See: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
        """

    def _strand_to_location_string(self) -> str:
        """
        Returns a string representation appropriate for location_string of the strand of the
        primer
        """
        match self.span.strand:
            case Strand.POSITIVE:
                return "F"
            case Strand.NEGATIVE:
                return "R"
            case _:  # pragma: no cover
                # Not calculating coverage on this line as it should be impossible to reach
                assert_never(f"Encountered unhandled Strand value: {self.span.strand}")


PrimerLikeType = TypeVar("PrimerLikeType", bound=PrimerLike)
"""Type variable for classes generic over `PrimerLike` types."""
