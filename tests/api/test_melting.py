import pytest

from prymer.api.melting import calculate_long_seq_tm


@pytest.mark.parametrize(
    "seq, salt_molar_concentration, percent_formamide, expected_tm",
    [
        ("A" * 100, 1.0, 0.0, 74.75),
        ("C" * 100, 1.0, 0.0, 115.75),
        ("G" * 100, 1.0, 0.0, 115.75),
        ("T" * 100, 1.0, 0.0, 74.75),
        ("AT" * 50, 1.0, 0.0, 74.75),
        ("GC" * 50, 1.0, 0.0, 115.75),
        ("AC" * 50, 1.0, 0.0, 95.25),
        ("GT" * 50, 1.0, 0.0, 95.25),
        ("GT" * 10, 1.0, 0.0, 68.25),
        ("GT" * 10, 1.0, 2.0, 67.01),
        ("GT" * 10, 1.0, 10.0, 62.05),
        ("GT" * 10, 2.0, 10.0, 67.0471),
        ("GT" * 10, 10.0, 10.0, 78.65),
    ],
)
def test_calculate_long_seq_tm(
    seq: str, salt_molar_concentration: float, percent_formamide: float, expected_tm: float
) -> None:
    assert pytest.approx(expected_tm) == calculate_long_seq_tm(
        seq=seq,
        salt_molar_concentration=salt_molar_concentration,
        percent_formamide=percent_formamide,
    )
