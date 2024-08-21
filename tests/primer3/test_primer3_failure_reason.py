from collections import Counter
from typing import Optional

import pytest

from prymer.primer3.primer3_failure_reason import Primer3FailureReason

_FAILURE_REASONS_AND_ENUM_VALUES: list[tuple[str, Primer3FailureReason]] = [
    ("GC content failed", Primer3FailureReason.GC_CONTENT),
    ("unknown failure reason", None),
    ("long dinucleotide run", Primer3FailureReason.LONG_DINUC),
    ("undesirable secondary structure", Primer3FailureReason.SECONDARY_STRUCTURE),
    ("amplifies off-target regions", Primer3FailureReason.OFF_TARGET_AMPLIFICATION),
    ("amplifies off-target regions", Primer3FailureReason.OFF_TARGET_AMPLIFICATION),
]


@pytest.fixture
def failure_counter() -> Counter:
    counter: Counter = Counter()
    count = 1
    for failure_reason, enum_value in _FAILURE_REASONS_AND_ENUM_VALUES:
        if enum_value is not None:  # for the "unknown failure reason" above
            counter[failure_reason] += count
        count += 1
    return counter


@pytest.mark.parametrize("failure_reason, expected_enum_value", _FAILURE_REASONS_AND_ENUM_VALUES)
def test_from_reason(
    failure_reason: str, expected_enum_value: Optional[Primer3FailureReason]
) -> None:
    assert Primer3FailureReason.from_reason(failure_reason) == expected_enum_value


@pytest.mark.parametrize("failure_reason, expected_enum_value", _FAILURE_REASONS_AND_ENUM_VALUES)
def test_parse_failures_individually(
    failure_reason: str, expected_enum_value: Optional[Primer3FailureReason]
) -> None:
    if expected_enum_value is None:
        # test that the failure reason cannot be parsed
        failure_string = failure_reason
        counter = Primer3FailureReason.parse_failures(failure_string)
        assert expected_enum_value not in counter
    else:
        failure_string = f"{failure_reason} {1}"
        counter = Primer3FailureReason.parse_failures(failure_string)
        assert counter[expected_enum_value] == 1


def test_parse_failures_combined(failure_counter: Counter) -> None:
    # build the failure strings
    failure_strings = []
    for failure_reason, count in failure_counter.items():
        if count % 2 == 0:
            # just have one reason with a count
            failure_string = f"{failure_reason} {count}"
            failure_strings.append(failure_string)
        else:
            # add the same reason _count_ # of times
            for _ in range(count):
                failure_string = f"{failure_reason} 1"
                failure_strings.append(failure_string)

    # test it giving one string at a time
    assert Primer3FailureReason.parse_failures(*failure_strings) == failure_counter

    # test it giving all the strings, comma-delimited
    assert Primer3FailureReason.parse_failures(", ".join(failure_strings)) == failure_counter
