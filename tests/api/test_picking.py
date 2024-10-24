import dataclasses
from typing import Optional

import pytest
from fgpyo.sequence import reverse_complement

from prymer.api import MinOptMax
from prymer.api import Oligo
from prymer.api import PrimerPair
from prymer.api import Span
from prymer.api import picking
from prymer.api.melting import calculate_long_seq_tm
from prymer.primer3 import PrimerAndAmpliconWeights


@pytest.fixture
def amplicon_sizes() -> MinOptMax[int]:
    return MinOptMax(100, 150, 200)


@pytest.fixture
def amplicon_tms() -> MinOptMax[float]:
    return MinOptMax(75.0, 80.0, 90.0)


@pytest.fixture
def weights() -> PrimerAndAmpliconWeights:
    return PrimerAndAmpliconWeights(
        product_size_lt=0.5,
        product_size_gt=1.5,
        product_tm_lt=10,
        product_tm_gt=15,
    )


@pytest.fixture
def all_zero_weights() -> PrimerAndAmpliconWeights:
    d = dataclasses.asdict(PrimerAndAmpliconWeights())
    for key in d:
        d[key] = 0.0
    return PrimerAndAmpliconWeights(**d)


def p(bases: str, tm: float, pos: int, pen: float = 0) -> Oligo:
    """Generates a primer for testing."""
    return Oligo(
        name="left",
        tm=tm,
        penalty=pen,
        bases=bases,
        span=Span("chr1", pos, pos + len(bases) - 1),
        tail=None,
    )


def pp(lp: Oligo, rp: Oligo, bases: Optional[str] = None, tm: Optional[float] = None) -> PrimerPair:
    """Generates a primer pair for testing."""
    if bases is None:
        length = rp.span.end - lp.span.start + 1
        needed = length - lp.span.length - rp.span.length
        bases = lp.bases + ("A" * needed) + reverse_complement(rp.bases)

    return PrimerPair(
        name="pair",
        left_primer=lp,
        right_primer=rp,
        amplicon_tm=tm if tm is not None else calculate_long_seq_tm(bases),
        amplicon_sequence=bases,
        penalty=0,
    )


def _score(
    pair: PrimerPair,
    weights: PrimerAndAmpliconWeights,
    sizes: MinOptMax[int],
    tms: MinOptMax[float],
) -> float:
    """Wrapper around picking.score that decomposes the PrimerPair for the arguments."""
    return picking.score(
        left_primer=pair.left_primer,
        right_primer=pair.right_primer,
        amplicon=pair.span,
        amplicon_tm=pair.amplicon_tm,
        amplicon_sizes=sizes,
        amplicon_tms=tms,
        weights=weights,
    )


def test_score_returns_sum_of_primer_penalties_when_all_weights_zero(
    all_zero_weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 100, 1.7), p("GTGCATACTAG", 59.8, 200, 3.1))
    score = _score(pair=pair, weights=all_zero_weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert score == 1.7 + 3.1


def test_score_returns_sum_of_primer_penalties_when_amplicon_optimal(
        all_zero_weights: PrimerAndAmpliconWeights,
        amplicon_sizes: MinOptMax[int],
        amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 100, 1.7), p("GTGCATACTAG", 59.8, 240, 3.1), tm=80)
    score = _score(pair=pair, weights=all_zero_weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert score == 1.7 + 3.1
