"""
# Utility Classes and Methods for NtThermoAlign

This module contains the [`NtThermoAlign`][prymer.ntthal.NtThermoAlign] class for
submitting multiple queries to `ntthal`, a command line tool included with `primer3`. Methods
include functions to open and use subprocesses, check and
manipulate status, and calculate melting temperature of valid oligo sequences.
The class can be optionally used as a context manager.

## Examples of Calculating Melting Temperature

```python
>>> from prymer import ntthal
>>> t = ntthal.NtThermoAlign()
>>> print(t.duplex_tm(s1="ATGC", s2="GCAT"))
-54.75042
>>> from prymer import ntthal
>>> with ntthal.NtThermoAlign() as t:
...     print(t.duplex_tm(s1="ATGC", s2="GCAT"))
-54.75042

```
"""

from pathlib import Path

from prymer.util.executable_runner import ExecutableRunner

MONOVALENT_MILLIMOLAR: float = 50.0
"""The default concentration of monovalent cations in mM"""

DIVALENT_MILLIMOLAR: float = 0.0
"""The default concentration of divalent cations in mM"""

DNTP_MILLIMOLAR: float = 0.0
"""The default concentration of deoxynycleotide triphosphate in mM"""

DNA_NANOMOLAR: float = 50.0
"""The concentration of DNA strands in nM"""

TEMPERATURE: float = 37.0
"""The default temperature at which duplex is calculated (Celsius)"""


class NtThermoAlign(ExecutableRunner):
    """
    Uses the `ntthal` command line tool to calculate the melting temperature of user-provided
    oligo sequences.
    """

    def __init__(
        self,
        executable: str | Path = "ntthal",
        monovalent_millimolar: float = MONOVALENT_MILLIMOLAR,
        divalent_millimolar: float = DIVALENT_MILLIMOLAR,
        dntp_millimolar: float = DNTP_MILLIMOLAR,
        dna_nanomolar: float = DNA_NANOMOLAR,
        temperature: float = TEMPERATURE,
    ):
        """
        Args:
            executable: string or Path representation of ntthal executable path
            monovalent_millimolar: concentration of monovalent cations in mM
            divalent_millimolar: concentration of divalent cations in mM
            dntp_millimolar: concentration of deoxynycleotide triphosphate in mM
            dna_nanomolar: concentration of DNA strands in nM
            temperature: temperature at which duplex is calculated (Celsius)
        """
        executable_path = ExecutableRunner.validate_executable_path(executable=executable)
        command: list[str] = [f"{executable_path}", "-r", "-i"]

        if monovalent_millimolar < 0:
            raise ValueError(f"monovalent_millimolar must be >=0, received {monovalent_millimolar}")
        if divalent_millimolar < 0:
            raise ValueError(f"divalent_millimolar must be >=0, received {divalent_millimolar}")
        if dntp_millimolar < 0:
            raise ValueError(f"dntp_millimolar must be >=0, received {dntp_millimolar}")
        if dna_nanomolar < 0:
            raise ValueError(f"dna_nanomolar must be >=0, received {dna_nanomolar}")
        if temperature < 0:
            raise ValueError(f"temperature must be >=0, received {temperature}")

        command.extend(["-mv", f"{monovalent_millimolar}"])
        command.extend(["-dv", f"{divalent_millimolar}"])
        command.extend(["-n", f"{dntp_millimolar}"])
        command.extend(["-d", f"{dna_nanomolar}"])
        command.extend(["-t", f"{temperature}"])

        super().__init__(command=command)

    def duplex_tm(self, s1: str, s2: str) -> float:
        """
        Calculates the melting temperature (Tm) of two provided oligos.

        Args:
            s1: the sequence of oligo 1 (5'->3' orientation)
            s2: the sequence of oligo 2 (5'->3' orientation)

        Example:
            >>> t = NtThermoAlign()
            >>> t.duplex_tm(s1 = "ACGT", s2 = "ACGT")
            -46.542706

        Returns:
            result: ntthal-calculated melting temperature

        Raises:
            ValueError: if ntthal result cannot be cast to a float
            RuntimeError: if underlying subprocess has already been terminated
        """
        if not self.is_alive:
            raise RuntimeError(
                "Error, trying to use a subprocess that has already been "
                f"terminated, return code {self._subprocess.returncode}"
            )
        if not s1.isalpha() or not s2.isalpha():
            raise ValueError(
                "Both input strings must be all alphabetic and non-empty, "
                f"received {s1} and {s2}"
            )

        self._subprocess.stdin.write(f"{s1},{s2}\n")
        self._subprocess.stdin.flush()  # forces the input to be sent to the underlying process.
        raw_result = self._subprocess.stdout.readline().rstrip("\r\n")
        try:
            result = float(raw_result)
        except ValueError as e:
            raise ValueError(f"Error: {e}, {raw_result} cannot be cast to float") from e

        return result
