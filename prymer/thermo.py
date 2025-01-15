from dataclasses import dataclass
from dataclasses import field
from typing import ClassVar

from primer3.thermoanalysis import ThermoResult  # type: ignore
from primer3.thermoanalysis import _ThermoAnalysis


@dataclass(frozen=True, kw_only=True)
class Thermo:
    """
    Class for performing thermodynamic calculations.  Largely wraps the `_ThermoAnalysis`
    class from `primer3-py` to provide a more ergonomic interface and reset default thermodynamic
    parameters to match those used in primer3.

    mv_conc_mm: concentration of monovalent cations in mM
    dv_conc_mm: concentration of divalent cations in mM
    dntp_conc_mm: concentration of deoxynycleotide triphosphate in mM
    dna_conc_nm: concentration of DNA strands in nM
    dmso_conc_pct: Concentration of DMSO (%)
    dmso_correction_factor: DMSO correction factor
    formamide_conc_mol: Concentration of formamide (mol/l)
    annealing_temp_c: Actual annealing temperature of the PCR reaction in Celsius
    temp_c: Simulation temperature for structure prediction in Celsius
    max_nn_length: Maximum length in bases for nearest-neighbor calcs
    max_loop: maximum size of loops in bases when predicting structures
    tm_method: Tm calculation method
    salt_corrections_method: Salt correction method
    """

    TM_METHOD_BRESLAUER: ClassVar[str] = "breslauer"
    TM_METHOD_SANTALUCIA: ClassVar[str] = "santalucia"

    SALT_CORRECTION_METHOD_SCHILDKRAUT: ClassVar[str] = "schildkraut"
    SALT_CORRECTION_METHOD_OWCZARZY: ClassVar[str] = "owczarzy"
    SALT_CORRECTION_METHOD_SANTALUCIA: ClassVar[str] = "santalucia"

    mv_conc_mm: float = 50.0
    dv_conc_mm: float = 0.0
    dntp_conc_mm: float = 0.0
    dna_conc_nm: float = 50.0
    dmso_conc_pct: float = 0.0
    dmso_correction_factor: float = 0.6
    formamide_conc_mol: float = 0.0
    annealing_temp_c: float = -10.0
    temp_c: float = 37.0
    max_nn_length: int = 60
    max_loop: int = 30
    tm_method: str = TM_METHOD_SANTALUCIA
    salt_corrections_method: str = SALT_CORRECTION_METHOD_SANTALUCIA

    _thermo: _ThermoAnalysis = field(init=False, default_factory=lambda: _ThermoAnalysis())

    def __post_init__(self) -> None:
        self._thermo.set_thermo_args(
            mv_conc=self.mv_conc_mm,
            dv_conc=self.dv_conc_mm,
            dntp_conc=self.dntp_conc_mm,
            dna_conc=self.dna_conc_nm,
            dmso_conc=self.dmso_conc_pct,
            dmso_fact=self.dmso_correction_factor,
            formamide_conc=self.formamide_conc_mol,
            annealing_temp_c=self.annealing_temp_c,
            temp_c=self.temp_c,
            max_nn_length=self.max_nn_length,
            max_loop=self.max_loop,
            tm_method=self.tm_method,
            salt_corrections_method=self.salt_corrections_method,
        )

    def tm(self, bases: str) -> float:
        return float(self._thermo.calc_tm(bases))

    def hairpin_tm(self, bases: str) -> float:
        result: ThermoResult = self._thermo.calc_hairpin(bases)
        return float(result.tm)

    def homodimer_tm(self, bases: str) -> float:
        result: ThermoResult = self._thermo.calc_homodimer(bases)
        return float(result.tm)

    def heterodimer_tm(self, bases1: str, bases2: str) -> float:
        result: ThermoResult = self._thermo.calc_heterodimer(bases1, bases2)
        return float(result.tm)

    def heterodimer_3p_tm(self, bases1: str, bases2: str) -> float:
        result: ThermoResult = self._thermo.calc_end_stability(bases1, bases2)
        return float(result.tm)
