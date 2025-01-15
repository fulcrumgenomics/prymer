from dataclasses import replace

from prymer import Thermo

# This module provides some very basic tests for the Thermo class on the assumption that
# the underlying functionality from ntthal is robust.


def test_all_tm_and_correction_methods() -> None:
    salt_methods = [
        Thermo.SALT_CORRECTION_METHOD_SANTALUCIA,
        Thermo.SALT_CORRECTION_METHOD_SCHILDKRAUT,
        Thermo.SALT_CORRECTION_METHOD_OWCZARZY,
    ]

    tm_methods = [Thermo.TM_METHOD_SANTALUCIA, Thermo.TM_METHOD_BRESLAUER]

    bases = "CGCTAATCGGGCTAGCTAG"

    for tm_method in tm_methods:
        for salt_method in salt_methods:
            thermo = Thermo(tm_method=tm_method, salt_corrections_method=salt_method)
            tm = thermo.tm(bases)
            assert tm > 0


def test_tm_sequence_length() -> None:
    thermo = Thermo()
    assert thermo.tm("ACGTT" * 4) > 0  # short, 20bp sequence
    assert thermo.tm("ACGTT" * 40) > 0  # long, 200bp sequence


def test_hairpin_tm() -> None:
    thermo = Thermo()
    tm = thermo.hairpin_tm("CCCTTTTTAAAGGG")
    assert 50 < tm < 70


def test_homodimer_tm() -> None:
    thermo = Thermo()
    tm1 = thermo.homodimer_tm("AAAAAAAAAAAAAAAAAAAA")  # no homodimer possibly
    tm2 = thermo.homodimer_tm("AAAAAAAAAATTTTTTTTTT")  # makes a perfect homodimer

    assert tm1 == 0
    assert tm2 > 0


def test_heterodimer_tm() -> None:
    thermo = Thermo()
    tm1 = thermo.heterodimer_tm("AAAAAAAAAA", "CCCCCCCCCC")
    tm2 = thermo.heterodimer_tm("AAAAAAAAAA", "TTTTTTTTTT")

    assert tm1 == 0
    assert tm2 > 0


def test_heterodimer_3p_anchored_tm() -> None:
    thermo = Thermo()
    tm1 = thermo. heterodimer_3p_anchored_tm("AAAAAAAAAA", "CCCCCCCCCC")
    tm2 = thermo. heterodimer_3p_anchored_tm("AAAAAAAAAA", "TTTTTTTTTT")

    assert tm1 == 0
    assert tm2 > 0


def test_dataclass_replace_works() -> None:
    thermo = Thermo()
    bases = "CGCTAATCGGGCTAGCTAG"
    tm1 = thermo.tm(bases)
    thermo = replace(thermo, dmso_conc_pct=10)  # adding DMSO should drop Tms
    tm2 = thermo.tm(bases)

    assert tm2 < tm1
