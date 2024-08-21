from dataclasses import dataclass
from typing import Optional

import pytest

from prymer.api.primer_like import PrimerLike
from prymer.api.span import Span


@dataclass(frozen=True, init=True, kw_only=True, slots=True)
class PrimerLikeTester(PrimerLike):
    """A simple class that inherits from PrimerLike for testing purposes."""

    span: Span
    bases: Optional[str] = None

    def to_bed12_row(self) -> str:
        return "a/b/c/d/e/f/g/h/i/j/k/l"


@pytest.mark.parametrize(
    "name, expected_id",
    [
        ("test", "test"),
        (None, "chr1_1_10_F"),
    ],
)
def test_id_generation(
    test_span: Span,
    name: Optional[str],
    expected_id: str,
) -> None:
    """Asserts that the id field is correctly generated based on the name and name_prefix fields."""
    test_primer = PrimerLikeTester(name=name, bases="AATCGATCCA", span=test_span)
    assert test_primer.id == expected_id


def test_to_bed12_row_exists(test_span: Span) -> None:
    """Asserts that the to_bed12_row method exists and returns the expected value."""
    test_primer = PrimerLikeTester(
        name="test",
        bases="AATCGATCCA",
        span=test_span,
    )
    assert test_primer.to_bed12_row() == "a/b/c/d/e/f/g/h/i/j/k/l"
