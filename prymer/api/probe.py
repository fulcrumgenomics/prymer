from dataclasses import dataclass

from fgpyo.util.metric import Metric

from prymer.api.primer import Primer


@dataclass(frozen=True, init=True, kw_only=True, slots=True)
class Probe(Primer, Metric["Probe"]):
    """Stores the properties of the designed Probe. Inherits `tm`, `penalty`,
    `span`, `bases`, and `tail` from `Primer`.

    Attributes:
        self_any_th: self-complementarity throughout the probe as calculated by Primer3
        self_end_th: 3' end complementarity of the probe as calculated by Primer3
        hairpin_th: hairpin formation thermodynamics of the probe as calculated by Primer3

    """

    self_any_th: float
    self_end_th: float
    hairpin_th: float
