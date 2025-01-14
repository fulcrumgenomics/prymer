"""
# Primer3FailureReason Class

This module contains a Primer3FailureReason class. Based on user-specified criteria, Primer3 will
disqualify primer designs if the characteristics of the design are outside an allowable range of
parameters.

Failure reason strings are documented in the Primer3 source code, accessible here
(at the time of Primer3FailureReason implementation):
https://github.com/bioinfo-ut/primer3_masker/blob/master/src/libprimer3.c#L5581-L5604

## Example Primer3 failure reasons emitted

Primer3 shall return a comma-delimited list of failure explanations.  They will have key
`PRIMER_LEFT_EXPLAIN` for designing individual left primers, `PRIMER_RIGHT_EXPLAIN` for designing
individual right primers, and `PRIMER_PAIR_EXPLAIN` for designing primer pairs.

For individual primers:

```python
>>> failure_string = 'considered 160, too many Ns 20, low tm 127, ok 13'
>>> Primer3FailureReason.parse_failures(failure_string)
Counter({<Primer3FailureReason.LOW_TM: 'low tm'>: 127, <Primer3FailureReason.TOO_MANY_NS: 'too many Ns'>: 20})
>>> failure_string = 'considered 238, low tm 164, high tm 12, high hairpin stability 23, ok 39'
>>> Primer3FailureReason.parse_failures(failure_string)
Counter({<Primer3FailureReason.LOW_TM: 'low tm'>: 164, <Primer3FailureReason.HAIRPIN_STABILITY: 'high hairpin stability'>: 23, <Primer3FailureReason.HIGH_TM: 'high tm'>: 12})
>>> failure_string = 'considered 166, unacceptable product size 161, ok 5'
>>> Primer3FailureReason.parse_failures(failure_string)
Counter({<Primer3FailureReason.PRODUCT_SIZE: 'unacceptable product size'>: 161})

```

"""  # noqa: E501

import logging
import re
from collections import Counter
from enum import StrEnum
from enum import unique
from typing import Optional


@unique
class Primer3FailureReason(StrEnum):
    """
    Enum to represent the various reasons Primer3 removes primers and primer pairs.

    These are taken from: https://github.com/bioinfo-ut/primer3_masker/blob/
      6a40c4c408dc02b95ac02391457cda760092291a/src/libprimer3.c#L5581-L5605

    This also contains custom failure values that are not generated by Primer3 but are convenient
    to have so that post-processing failures can be tracked in the same way that Primer3 failures
    are.  These include `LONG_DINUC`, `SECONDARY_STRUCTURE`, and `OFF_TARGET_AMPLIFICATION`.
    """

    # Failure reasons emitted by Primer3
    GC_CONTENT = "GC content failed"
    GC_CLAMP = "GC clamp failed"
    HAIRPIN_STABILITY = "high hairpin stability"
    HIGH_TM = "high tm"
    LOW_TM = "low tm"
    LOWERCASE_MASKING = "lowercase masking of 3' end"
    LONG_POLY_X = "long poly-x seq"
    PRODUCT_SIZE = "unacceptable product size"
    TOO_MANY_NS = "too many Ns"
    HIGH_ANY_COMPLEMENTARITY = "high any compl"
    HIGH_END_COMPLEMENTARITY = "high end compl"
    # Failure reasons not emitted by Primer3 and beneficial to keep track of
    LONG_DINUC = "long dinucleotide run"
    SECONDARY_STRUCTURE = "undesirable secondary structure"
    OFF_TARGET_AMPLIFICATION = "amplifies off-target regions"

    @classmethod
    def from_reason(cls, str_reason: str) -> Optional["Primer3FailureReason"]:
        """Returns the first `Primer3FailureReason` with the given reason for failure.
        If no failure exists, return `None`."""
        reason: Optional[Primer3FailureReason] = None
        try:
            reason = cls(str_reason)
        except ValueError:
            pass
        return reason

    @staticmethod
    def parse_failures(
        *failures: str,
    ) -> Counter["Primer3FailureReason"]:
        """When Primer3 encounters failures, extracts the reasons why designs that were considered
         by Primer3 failed.

        Args:
            failures: list of strings, with each string an "explanation" emitted by Primer3 about
                why the design failed.  Each string may be a comma delimited of failures, or a
                single failure.

        Returns:
            a `Counter` of each `Primer3FailureReason` reason.
        """

        failure_regex = r"^ ?(.+) ([0-9]+)$"
        by_fail_count: Counter[Primer3FailureReason] = Counter()
        # parse all failure strings and merge together counts for the same kinds of failures
        for failure in failures:
            split_failures = [fail.strip() for fail in failure.split(",")]
            for item in split_failures:
                result = re.match(failure_regex, item)
                if result is None:
                    continue
                reason = result.group(1)
                count = int(result.group(2))
                if reason in ["ok", "considered"]:
                    continue
                std_reason = Primer3FailureReason.from_reason(reason)
                if std_reason is None:
                    logging.getLogger(__name__).debug(f"Unknown Primer3 failure reason: {reason}")
                by_fail_count[std_reason] += count

        return by_fail_count
