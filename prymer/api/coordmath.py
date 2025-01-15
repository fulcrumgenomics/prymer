"""
# Methods for coordinate-based math.
"""


def get_closed_end(start: int, length: int) -> int:
    """
    Get the closed end of an interval given its start and length.

    Args:
        start: The start position of the interval.
        length: The length of the interval.

    Returns:
        The closed end of the interval.

    Example:
        ```python
        >>> get_closed_end(start=10, length=5)
        14

        ```
    """
    return start + length - 1
