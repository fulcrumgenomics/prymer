from dataclasses import dataclass
from dataclasses import field
from typing import ClassVar

from primer3.thermoanalysis import ThermoResult  # type: ignore
from primer3.thermoanalysis import _ThermoAnalysis


@dataclass(frozen=True, kw_only=True)
class Thermo:
    """
    Class for performing thermodynamic calculations.  Available calculations include:

      1. melting temperature (Tm) for short and long sequences ([`tm()`][prymer.thermo.Thermo.tm])
      2. hairpin / secondary structure Tm for single sequences ([`hairpin_tm`][prymer.thermo.Thermo.hairpin_tm])
      3. homodimer Tm - the Tm of a duplex formed from two copies of the same sequence ([`homodimer_tm `][prymer.thermo.Thermo. homodimer_tm])
      4. heterodimer Tm - the Tm of a duplex formed from two different sequences ([`heterodimer_tm `][prymer.thermo.Thermo. heterodimer_tm])
      5. 3' anchored heterodimer Tm - the heterodimer Tm when annealing of the 3' end is prioritized ([`heterodimer_3p_anchored_tm `][prymer.thermo.Thermo. heterodimer_3p_anchored_tm])

    The `tm` method can be used for sequences of any length.  For sequences up to `max_nn_length`
    (default value 60 bases) the nearest neighborhood Tm calculation is used.  For sequences above
    the `max_nn_length` a more generic long sequence Tm calculation is used.

    Many attributes affect the calculation of Tms and these are set to the default values used by
    primer3 and, where applicable, ntthal.

    Class level constants are provided for Tm calculation method names and salt correction method
    names, and it is strongly suggested to use these.

    Attributes:
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
    dv_conc_mm: float = 1.5
    dntp_conc_mm: float = 0.6
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
        """Transfers the values to the underlying implementation class."""
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
        """
        Calculates the melting temperature of a sequence of bases.  The calculation used differs
        for short sequences (len < `max_nn_length`) and long sequences.

        Arguments:
            bases: a sequence of unambiguous DNA bases (i.e. A/C/G/T)

        Returns:
            the melting temperature of the sequence, in degrees Celsius
        """
        return float(self._thermo.calc_tm(bases))

    def hairpin_tm(self, bases: str) -> float:
        """
        Calculates the melting temperature of the most likely secondary structure or hairpin of a
        single oligo.

        Arguments:
            bases: a sequence of unambiguous DNA bases (i.e. A/C/G/T)

        Returns:
            the melting temperature of the most likely (lowest energy) secondary structure
        """
        result: ThermoResult = self._thermo.calc_hairpin(bases)
        return float(result.tm)

    def homodimer_tm(self, bases: str) -> float:
        """
        Calculates the melting temperature of the most likely annealing of a single oligo to a
        second copy of the same oligo (i.e. a duplex of two of the same oligo).

        Arguments:
            bases: a sequence of unambiguous DNA bases (i.e. A/C/G/T)

        Returns:
            the melting temperature of the most likely (lowest energy) duplex
        """
        result: ThermoResult = self._thermo.calc_homodimer(bases)
        return float(result.tm)

    def heterodimer_tm(self, bases1: str, bases2: str) -> float:
        """
        Calculates the melting temperature of the most likely annealing of two oligos.  The two
        oligos do *not* need to be the same length.

        Arguments:
            bases1: a sequence of unambiguous DNA bases (i.e. A/C/G/T) for the first oligo
            bases2: a sequence of unambiguous DNA bases (i.e. A/C/G/T) for the second oligo

        Returns:
            the melting temperature of the most likely (lowest energy) heterodimer
        """

        result: ThermoResult = self._thermo.calc_heterodimer(bases1, bases2)
        return float(result.tm)

    def heterodimer_3p_anchored_tm(self, bases1: str, bases2: str) -> float:
        """
        Calculates the melting temperature of the most likely annealing of two oligos where the
        annealing is constrained to optimize annealing of the 3' end of the first sequence.

        Arguments:
            bases1: a sequence of unambiguous DNA bases (i.e. A/C/G/T) for the first oligo
            bases2: a sequence of unambiguous DNA bases (i.e. A/C/G/T) for the second oligo

        Returns:
            the melting temperature of the most likely (lowest energy) 3' anchored heterodimer
        """
        result: ThermoResult = self._thermo.calc_end_stability(bases1, bases2)
        return float(result.tm)
