"""This module is deprecated - see prymer/api/oligo.py"""

import warnings
from dataclasses import dataclass

from prymer.api.oligo import Oligo


@dataclass(frozen=True, init=True, slots=True)
class Primer(Oligo):
    """
    A deprecated alias for `Oligo`.

    This class exists to maintain backwards compatibility with earlier releases of `prymer`
    and may be removed in a future version.
    """

    warnings.warn(
        "The Primer class was deprecated, use Oligo instead",
        DeprecationWarning,
        stacklevel=2,
    )
