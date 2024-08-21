"""
# Methods for calculating melting temperatures.


There is currently one public method available:

- [`calculate_long_seq_tm()`][prymer.api.melting.calculate_long_seq_tm] -- Calculates the
    melting temperature of an amplicon.

## Examples


```python
>>> calculate_long_seq_tm(seq="GT" * 10, salt_molar_concentration=10.0, percent_formamide=10.0)
78.64999999999999

```
"""

import math

from fgpyo.sequence import gc_content


def calculate_long_seq_tm(
    seq: str, salt_molar_concentration: float = 1.65, percent_formamide: float = 15.0
) -> float:
    """Calculate the melting temperature of an amplicon.

    Uses the formula:

    `Tm = 81.5 + 0.41(%GC) - 675/N + 16.6 x log[Na+] - 0.62(%F)`

    from:

    (Marmur & Doty 1962, J Mol Biol 5: 109-118; Schildkraut & Lifson 1965, Biopolymers 3: 195-208)

    with the added chemical (formamide) correction.

    Args:
        seq: the amplicon sequence
        salt_molar_concentration: the molar concentration of salt
        percent_formamide: the percent formamide

    Returns:
        the predicted melting temperature
    """
    return (
        81.5
        + (16.6 * math.log10(salt_molar_concentration))
        + (41.0 * gc_content(seq))
        - (675.0 / len(seq))
        - (0.62 * percent_formamide)
    )
