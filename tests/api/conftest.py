import pytest
from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.fasta.sequence_dictionary import SequenceMetadata


@pytest.fixture
def seq_dict() -> SequenceDictionary:
    metadatas: list[SequenceMetadata] = [
        SequenceMetadata(name="chr1", length=1000000, index=0),
        SequenceMetadata(name="chr2", length=1000000, index=1),
        SequenceMetadata(name="chr3", length=1000000, index=2),
    ]
    return SequenceDictionary(metadatas)
