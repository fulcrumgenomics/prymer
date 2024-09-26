"""
# Primer Class and Methods

This module contains a class and class methods to represent a primer (e.g. designed by Primer3)

Class attributes include the primer sequence, melting temperature, and the score of the primer. The
mapping of the primer to the genome is also stored.

Optional attributes include naming information and a tail sequence  to attach to the 5' end of the
primer (if applicable).

## Examples of interacting with the `Primer` class

```python
>>> from prymer.api.span import Span, Strand
>>> primer_span = Span(refname="chr1", start=1, end=20)
>>> primer = Oligo(tm=70.0, penalty=-123.0, span=primer_span)
>>> primer.longest_hp_length()
0
>>> primer.length
20
>>> primer.name is None
True
>>> primer = Oligo(tm=70.0, penalty=-123.0, span=primer_span, bases="GACGG"*4)
>>> primer.longest_hp_length()
3
>>> primer.untailed_length()
20
>>> primer.tailed_length()
20
>>> primer = primer.with_tail(tail="GATTACA")
>>> primer.untailed_length()
20
>>> primer.tailed_length()
27
>>> primer = primer.with_name(name="foobar")
>>> primer.name
'foobar'

```

Primers may also be written to a file and subsequently read back in, as the `Primer` class is an
`fgpyo` `Metric` class:

```python
>>> from pathlib import Path
>>> left_span = Span(refname="chr1", start=1, end=20)
>>> left = Oligo(tm=70.0, penalty=-123.0, span=left_span, bases="G"*20)
>>> right_span = Span(refname="chr1", start=101, end=120)
>>> right = Oligo(tm=70.0, penalty=-123.0, span=right_span, bases="T"*20)
>>> path = Path("/tmp/path/to/primers.txt")
>>> Oligo.write(path, left, right)  # doctest: +SKIP
>>> primers = Oligo.read(path)  # doctest: +SKIP
>>> list(primers)  # doctest: +SKIP
[
    Primer(tm=70.0, penalty=-123.0, span=amplicon_span, bases="G"*20),
    Primer(tm=70.0, penalty=-123.0, span=amplicon_span, bases="T"*20)
]

```
"""

import warnings
from dataclasses import dataclass

from prymer.api.oligo import Oligo


@dataclass(frozen=True, init=True, slots=True)
class Primer(Oligo):
    """A deprecated alias for `Oligo` intended to maintain backwards
    compatibility with earlier releases of `prymer`."""

    warnings.warn(
        "The Primer class was deprecated, use Oligo instead",
        DeprecationWarning,
        stacklevel=2,
    )
