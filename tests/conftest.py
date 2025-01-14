"""
Fixtures intended to be shared across multiple files in the tests directory.
"""

from pathlib import Path

import pytest
from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.fasta.sequence_dictionary import SequenceMetadata

from prymer.model import Span
from prymer.model import Strand


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


@pytest.fixture
def seq_dict() -> SequenceDictionary:
    metadatas: list[SequenceMetadata] = [
        SequenceMetadata(name="chr1", length=1000000, index=0),
        SequenceMetadata(name="chr2", length=1000000, index=1),
        SequenceMetadata(name="chr3", length=1000000, index=2),
    ]
    return SequenceDictionary(metadatas)
