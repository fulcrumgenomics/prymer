"""
# Methods for coordinate-based math and interval manipulation.

Contains the following public methods:

- [`require_same_refname()`][prymer.api.coordmath.require_same_refname] -- ensures that all
    provided intervals have the same reference name.
- [`get_locus_string()`][prymer.api.coordmath.get_locus_string] -- returns a formatted
    string for an interval (`<refname>:<start>-<end>`).
- [`get_closed_end()`][prymer.api.coordmath.get_closed_end] -- gets the closed end of an
    interval given its start and length.

"""

from pybedlite.overlap_detector import Interval


def require_same_refname(*intervals: Interval) -> None:
    """
    Require that the input intervals all have the same refname.

    Args:
        intervals: one or more intervals

    Raises:
        ValueError: if the intervals do not all have the same refname.
    """

    refnames = set(i.refname for i in intervals)
    if len(refnames) != 1:
        raise ValueError(f"`intervals` must have exactly one refname\n Found {sorted(refnames)}")


def get_locus_string(record: Interval) -> str:
    """
    Get the locus-string for an interval.

    The output string will have the format `<refname>:<start>-<end>`
    No conversion on coordinates is performed, so the output is 0-based, open-ended.

    Args:
        record: The interval to get the locus-string for.

    Returns:
        A locus-string for the interval.
    """
    return f"{record.refname}:{record.start}-{record.end}"


def get_closed_end(start: int, length: int) -> int:
    """
    Get the closed end of an interval given its start and length.

    Args:
        start: The start position of the interval.
        length: The length of the interval.

    Returns:
        The closed end of the interval.

    Example:
        >>> get_closed_end(start=10, length=5)
        14
    """
    return start + length - 1
