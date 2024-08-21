import pytest

from prymer import ntthal


@pytest.mark.parametrize(
    "s1,s2,tm",
    [
        ("AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC", 0.0),  # no match
        (
            "AAAACCCGCTTTGCTAGCTACG",
            "AACCCGCTTGCTAGCAAAAAAAAAACCCGCTTGCTAGCAAAAAAAAAACCCGCTTGCTAGC"
            "AAAAAAAAAACCCGCTTGCTAGCAAAAAAAAAACCCGCTTGCTAGCAA",
            23.974756,  # significantly different lengths
        ),
        ("AAAACCCGCTTTGCTAGCTACG", "CGTAGCTAGCAAAGCGGGTTTT", 55.903276),  # reverse complements
        ("AAAACCCGCTTTGCTAGCTACGNNNNN", "NNNNNCGTAGCTAGCAAAGCGGGTTTT", 55.903276),  # NNNs
        ("GTCAGTCA", "ACGTACGT", -89.739318),  # negative Tm return
        ("XYZZZZZZZZZNNNNNNN", "QUEEEEEEENNNNN", 0.000000),  # non-ACGT characters
        ("A", "T", -437.071890),  # 1 nt,
    ],
)
def test_duplex_tm_default_params(s1: str, s2: str, tm: float) -> None:
    """Test that NtThermoAlign().duplexTm() will return Tm across valid conditions"""
    with ntthal.NtThermoAlign() as t:
        assert t.duplex_tm(s1, s2) == pytest.approx(tm)


@pytest.mark.parametrize(
    "invalid_s1,invalid_s2",
    [
        ("A123", "T456"),
        ("", ""),
        ("--A", "--T"),
        ("AGCT", ""),
        ("", "ACGT"),
    ],
)
def test_duplex_tm_invalid_input_raises(invalid_s1: str, invalid_s2: str) -> None:
    """Test that NtThermoAlign().duplexTm() will raise an error with invalid inputs"""
    t = ntthal.NtThermoAlign()
    with pytest.raises(ValueError):
        t.duplex_tm(invalid_s1, invalid_s2)


def test_invalid_path() -> None:
    """Test that NtThermoAlign().duplexTm() will raise an error if given an invalid path"""
    with pytest.raises(ValueError):
        ntthal.NtThermoAlign(executable="invalid_path")


@pytest.mark.parametrize(
    "monovalent_millimolar,divalent_millimolar,dntp_millimolar,dna_nanomolar,temp,expected_tm",
    [
        (50.0, 0.0, 0.8, 50.0, 37.0, -54.750420),
        (500.0, 0.0, 0.8, 50.0, 37.0, -49.372162),
        (50.0, 2.0, 0.8, 50.0, 37.0, -51.771985),
        (50.0, 0.0, 2.0, 50.0, 37.0, -54.750420),
        (50.0, 0.0, 0.8, 2.0, 37.0, -67.205211),
        (50.0, 0.0, 0.8, 50.0, 500.0, -54.750420),
    ],
)
def test_valid_ntthal_params(
    monovalent_millimolar: float,
    divalent_millimolar: float,
    dntp_millimolar: float,
    dna_nanomolar: float,
    temp: float,
    expected_tm: float,
) -> None:
    """Test class instantiation with valid params and context manager methods"""
    with ntthal.NtThermoAlign(
        monovalent_millimolar=monovalent_millimolar,
        divalent_millimolar=divalent_millimolar,
        dntp_millimolar=dntp_millimolar,
        dna_nanomolar=dna_nanomolar,
        temperature=temp,
    ) as t:
        assert t.duplex_tm("ATGC", "GCAT") == pytest.approx(expected_tm)


def test_ntthal_teardown() -> None:
    """Test teardown of NtThermoAlign() with context manager methods"""
    with ntthal.NtThermoAlign() as t:
        assert t.is_alive
    assert not t.is_alive
    with pytest.raises(RuntimeError):
        t.duplex_tm(s1="ACGT", s2="ACGT")


@pytest.mark.parametrize(
    "monovalent_millimolar,divalent_millimolar,dntp_millimolar,dna_nanomolar,temp",
    [
        (-50.0, 0.0, 0.8, 0.0, 37.0),
        (0.0, -50.0, 0.8, 0.0, 37.0),
        (0.0, 0.0, -50.0, 0.0, 37.0),
        (0.0, 0.0, 0.8, -50.0, 37.0),
        (0.0, 0.0, 0.8, 0.0, -50.0),
    ],
)
def test_invalid_ntthal_params(
    monovalent_millimolar: float,
    divalent_millimolar: float,
    dntp_millimolar: float,
    dna_nanomolar: float,
    temp: float,
) -> None:
    with pytest.raises(ValueError):
        ntthal.NtThermoAlign(
            monovalent_millimolar=monovalent_millimolar,
            divalent_millimolar=divalent_millimolar,
            dntp_millimolar=dntp_millimolar,
            dna_nanomolar=dna_nanomolar,
            temperature=temp,
        )
