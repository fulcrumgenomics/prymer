from pathlib import Path

import pytest


@pytest.fixture
def ref_fasta() -> Path:
    return Path(__file__).parent / "data" / "miniref.fa"
