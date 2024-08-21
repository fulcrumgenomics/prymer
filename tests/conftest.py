"""
Fixtures intended to be shared across multiple files in the tests directory.
"""

from pathlib import Path

import pytest

from prymer.api.span import Span
from prymer.api.span import Strand


@pytest.fixture
def test_span() -> Span:
    """A basic Span for use in tests that don't depend on specific values"""
    return Span(refname="chr1", start=1, end=10, strand=Strand.POSITIVE)


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def genome_ref(data_dir: Path) -> Path:
    return data_dir / "miniref.fa"
