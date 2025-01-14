import dataclasses
from math import floor
from pathlib import Path
from typing import Optional

import pytest
from fgpyo.fasta.builder import FastaBuilder
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


# 1000 bases with 100 bases per line
REF_BASES = """
CCTGCTGGTTTTATTTCTGCTACTTCTGAAGTGTGGGGCACACAACCTGAGCAGTGCTTTTATTTGAGTCCCAATGCTTTTATTTGAGTTTTGCAAGGTT
ATTCCAAGTTTTACAAATAGAAGGTAGCGTATGACTCAGTCCTTGATATGCCAACCACTGCACAGAGACTTGCCACCTTCCTGTCACTGGAGAAACACTC
ATGTGGGTTTTCTTAAATTTGCCTCCCTCTGAGCTTCCCTTTAACTTCAACTATAATATGCAAGAAAGACTATCTGACCATAAATACACATTTGGGCCAA
TCAAGATGGTTTTGCCAAGGAAAGATGCCCACAATGGTTAAGCAGAATGCAATAATGTAGAGAATATCATTTCTTTCATGCTGGTGTATATCATATGCAT
TCAAAAACAGGGAGAACTTCTAAGCAACTAACAGTGACCATATCAAGCAGGTGCAATCACAGAATAACTGGTTTTCTCCTTTAAGAATTTTTCTATCATT
TGGCTTTCCCCACTCACACACACTAAATATTTTAAGTAAAAAGTTACTTCCATTTTGAAAGAGAAAAGAAAGAGACATGCATGAACATTTTTCTCCACCT
TGGTGCAGGGACCAGACAACTGTATCCAGTGTGCCCACTACATTGACGGCCCCCACTGCGTCAAGACCTGCCCGGCAGGAGTCATGGGAGAAAACAACAC
CCTGGTCTGGAAGTACGCAGACGCCGGCCATGTGTGCCACCTGTGCCATCCAAACTGCACCTACGGGTGAGTGGAAAGTGAAGGAGAACAGAACATTTCC
TCTCTTGCAAATTCAGAGATCAAAAATGTCTCCCAAGTTTTCCGGCAACAAATTGCCGAGGTTTGTATTTGAGTCAGTTACTTAAGGTGTTTTGGTCCCC
ACAGCCATGCCAGTAGCAACTTGCTTGTGAGCAGGCCTCAGTGCAGTGGGAATGACTCTGCCATGCACCGTGTCCCCGGCCGGGCCTGTGTTGTGCAATG
""".strip().replace("\n", "")


@pytest.fixture
def fasta(tmp_path: Path) -> Path:
    """Fixture that returns a Fasta"""
    builder = FastaBuilder(line_length=100)
    path = tmp_path / "ref.fa"
    builder.add("chr1").add(REF_BASES)
    builder.to_file(path)
    return path


def p(bases: str, tm: float, pos: int, pen: float = 0, chrom: str = "chr1") -> Oligo:
    """Generates a primer for testing."""
    oligo = Oligo(
        name="left",
        tm=tm,
        penalty=pen,
        bases=bases,
        span=Span(chrom, pos, pos + len(bases) - 1),
        tail=None,
    )
    assert oligo.span.length == len(oligo.bases)
    return oligo


def pp(lp: Oligo, rp: Oligo, bases: Optional[str] = None, tm: Optional[float] = None) -> PrimerPair:
    """Generates a primer pair for testing."""
    if lp.span.end > rp.span.start:
        raise ValueError("Overlapping primers not supported.")

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


################################################################################
# Tests for picking.py::score()
################################################################################


def test_score_returns_sum_of_primer_penalties_when_all_weights_zero(
    all_zero_weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 100, 1.7), p("GTGCATACTAG", 59.8, 200, 3.1))
    score = _score(pair=pair, weights=all_zero_weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert score == 1.7 + 3.1


def test_score_returns_sum_of_primer_penalties_when_amplicon_optimal(
    weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 101, 1.7), p("GTGCATACTAG", 59.8, 240, 3.1), tm=80)
    score = _score(pair=pair, weights=weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert pair.amplicon.length == amplicon_sizes.opt
    assert pair.amplicon_tm == amplicon_tms.opt
    assert score == 1.7 + 3.1


def test_score_when_amplicon_longer_than_optimal(
    weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 101, 1.0), p("GTGCATACTAG", 59.8, 250, 1.0), tm=80)
    score = _score(pair=pair, weights=weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert pair.amplicon.length == amplicon_sizes.opt + 10
    assert pair.amplicon_tm == amplicon_tms.opt
    assert score == pytest.approx(1 + 1 + (10 * weights.product_size_gt))


def test_score_when_amplicon_shorter_than_optimal(
    weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 101, 1.0), p("GTGCATACTAG", 59.8, 220, 1.0), tm=80)
    score = _score(pair=pair, weights=weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert pair.amplicon.length == amplicon_sizes.opt - 20
    assert pair.amplicon_tm == amplicon_tms.opt
    assert score == pytest.approx(1 + 1 + (20 * weights.product_size_lt))


def test_score_when_amplicon_tm_higher_than_optimal(
    weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 101, 1.0), p("GTGCATACTAG", 59.8, 240, 1.0), tm=82.15)
    score = _score(pair=pair, weights=weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert pair.amplicon.length == amplicon_sizes.opt
    assert pair.amplicon_tm == amplicon_tms.opt + 2.15
    assert score == pytest.approx(1 + 1 + (2.15 * weights.product_tm_gt))


def test_score_when_amplicon_tm_lower_than_optimal(
    weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 101, 1.0), p("GTGCATACTAG", 59.8, 240, 1.0), tm=75.67)
    score = _score(pair=pair, weights=weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert pair.amplicon.length == amplicon_sizes.opt
    assert pair.amplicon_tm == amplicon_tms.opt - 4.33
    assert score == pytest.approx(1 + 1 + (4.33 * weights.product_tm_lt))


def test_score_realistic(
    weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pair = pp(p("ACGACTCATG", 60.1, 101, 3.1), p("GTGCATACTAG", 59.8, 256, 4.08), tm=75.67)
    score = _score(pair=pair, weights=weights, sizes=amplicon_sizes, tms=amplicon_tms)
    assert pair.amplicon.length == amplicon_sizes.opt + 16
    assert pair.amplicon_tm == amplicon_tms.opt - 4.33
    assert score == pytest.approx(
        3.1 + 4.08 + (16 * weights.product_size_gt) + (4.33 * weights.product_tm_lt)
    )


################################################################################
# Tests for picking.py::build_primer_pairs()
################################################################################


def test_build_primer_pairs_no_primers(
    fasta: Path,
    weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pairs = list(
        picking.build_primer_pairs(
            left_primers=[],
            right_primers=[],
            target=Span("chr1", 240, 260),
            amplicon_sizes=amplicon_sizes,
            amplicon_tms=amplicon_tms,
            max_heterodimer_tm=None,
            weights=weights,
            fasta_path=fasta,
        )
    )
    assert len(pairs) == 0


def test_build_primer_pairs_single_pair(
    fasta: Path,
    weights: PrimerAndAmpliconWeights,
    amplicon_sizes: MinOptMax[int],
    amplicon_tms: MinOptMax[float],
) -> None:
    pairs = list(
        picking.build_primer_pairs(
            left_primers=[p(REF_BASES[200:220], tm=60.1, pos=201, pen=1.0)],
            right_primers=[p(reverse_complement(REF_BASES[280:300]), tm=61.8, pos=281, pen=1.1)],
            target=Span("chr1", 240, 260),
            amplicon_sizes=amplicon_sizes,
            amplicon_tms=amplicon_tms,
            max_heterodimer_tm=None,
            weights=weights,
            fasta_path=fasta,
        )
    )

    for pair in pairs:
        assert pair.span.length == len(pair.bases)
        assert pair.amplicon.start == pair.left_primer.span.start
        assert pair.amplicon.end == pair.right_primer.span.end
        assert pair.bases == REF_BASES[pair.amplicon.start - 1 : pair.amplicon.end]

    assert len(pairs) == 1
    assert pairs[0].span == Span("chr1", 201, 300)
    assert pairs[0].left_primer.bases == pairs[0].bases[0:20]
    assert pairs[0].right_primer.bases == reverse_complement(pairs[0].bases[-20:])


def test_build_primers_amplicon_size_filtering(
    fasta: Path,
    weights: PrimerAndAmpliconWeights,
) -> None:
    pairs = list(
        picking.build_primer_pairs(
            left_primers=[p("A" * 20, tm=60, pos=s, pen=1.0) for s in range(51, 301, 10)],
            right_primers=[p("A" * 20, tm=60, pos=s, pen=1.0) for s in range(1, 401, 10)],
            target=Span("chr1", 150, 160),
            amplicon_sizes=MinOptMax(100, 150, 200),
            amplicon_tms=MinOptMax(0.0, 0.0, 100.0),
            max_heterodimer_tm=None,
            weights=weights,
            fasta_path=fasta,
        )
    )

    assert len(pairs) > 0

    for pair in pairs:
        assert pair.span.length == len(pair.bases)
        assert pair.amplicon.start == pair.left_primer.span.start
        assert pair.amplicon.end == pair.right_primer.span.end
        assert pair.bases == REF_BASES[pair.amplicon.start - 1 : pair.amplicon.end]
        assert pair.span.length <= 200
        assert pair.span.length >= 100


def test_build_primers_heterodimer_filtering(
    fasta: Path,
    weights: PrimerAndAmpliconWeights,
) -> None:
    pairs = list(
        picking.build_primer_pairs(
            left_primers=[
                p("CGCGCGCGCGCGCGCGCGCG", tm=60, pos=100, pen=1.0),
                p("CTGCTGCTGCTGCTGCTGCT", tm=60, pos=100, pen=1.0),
            ],
            right_primers=[
                p("GCGCGCGCGCGCGCGCGCGC", tm=60, pos=180, pen=1.0),
                p("AGCAGCAGCAGCAGCAGCAG", tm=60, pos=180, pen=1.0),
            ],
            target=Span("chr1", 150, 160),
            amplicon_sizes=MinOptMax(10, 150, 200),
            amplicon_tms=MinOptMax(0.0, 0.0, 100.0),
            max_heterodimer_tm=50,
            weights=weights,
            fasta_path=fasta,
        )
    )

    for pair in pairs:
        assert pair.span.length == len(pair.bases)
        assert pair.amplicon.start == pair.left_primer.span.start
        assert pair.amplicon.end == pair.right_primer.span.end
        assert pair.bases == REF_BASES[pair.amplicon.start - 1 : pair.amplicon.end]
        assert pair.span.length <= 200

    assert len(pairs) == 2
    primer_seqs = sorted([f"{pp.left_primer.bases}-{pp.right_primer.bases}" for pp in pairs])
    assert primer_seqs == [
        "CGCGCGCGCGCGCGCGCGCG-AGCAGCAGCAGCAGCAGCAG",
        "CTGCTGCTGCTGCTGCTGCT-GCGCGCGCGCGCGCGCGCGC",
    ]


def test_build_primer_pairs_amplicon_tm_filtering(
    fasta: Path,
    weights: PrimerAndAmpliconWeights,
) -> None:
    amp_bases = REF_BASES[200:300]
    amp_tm = calculate_long_seq_tm(amp_bases)
    assert floor(amp_tm) == 84

    for max_tm in [83, 84, 85]:
        pairs = list(
            picking.build_primer_pairs(
                left_primers=[p(REF_BASES[200:220], tm=60.1, pos=201, pen=1.0)],
                right_primers=[
                    p(reverse_complement(REF_BASES[280:300]), tm=61.8, pos=281, pen=1.1)
                ],
                target=Span("chr1", 240, 260),
                amplicon_sizes=MinOptMax(0, 0, 500),
                amplicon_tms=MinOptMax(0, 80, max_tm),
                max_heterodimer_tm=None,
                weights=weights,
                fasta_path=fasta,
            )
        )

        assert len(pairs) == (1 if max_tm > amp_tm else 0)


def test_build_primer_pairs_fails_when_primers_on_wrong_reference(
    fasta: Path,
    weights: PrimerAndAmpliconWeights,
) -> None:
    target = Span("chr1", 240, 260)
    valid_lefts = [p(REF_BASES[200:220], tm=60, pos=201, pen=1)]
    invalid_lefts = [p(REF_BASES[200:220], tm=60, chrom="X", pos=201, pen=1)]
    valid_rights = [p(reverse_complement(REF_BASES[280:300]), tm=61, pos=281, pen=1)]
    invalid_rights = [p(reverse_complement(REF_BASES[280:300]), tm=61, chrom="X", pos=281, pen=1)]

    picks = picking.build_primer_pairs(
        left_primers=valid_lefts,
        right_primers=valid_rights,
        target=target,
        amplicon_sizes=MinOptMax(0, 100, 500),
        amplicon_tms=MinOptMax(0, 80, 150),
        max_heterodimer_tm=None,
        weights=weights,
        fasta_path=fasta,
    )

    assert next(picks) is not None

    with pytest.raises(ValueError, match="Left primers exist on different reference"):
        _picks = list(
            picking.build_primer_pairs(
                left_primers=invalid_lefts,
                right_primers=valid_rights,
                target=target,
                amplicon_sizes=MinOptMax(0, 100, 500),
                amplicon_tms=MinOptMax(0, 80, 150),
                max_heterodimer_tm=None,
                weights=weights,
                fasta_path=fasta,
            )
        )

    with pytest.raises(ValueError, match="Right primers exist on different reference"):
        _picks = list(
            picking.build_primer_pairs(
                left_primers=valid_lefts,
                right_primers=invalid_rights,
                target=target,
                amplicon_sizes=MinOptMax(0, 100, 500),
                amplicon_tms=MinOptMax(0, 80, 150),
                max_heterodimer_tm=None,
                weights=weights,
                fasta_path=fasta,
            )
        )
