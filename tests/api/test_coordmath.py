import pytest
from pybedlite.overlap_detector import Interval

from prymer.api.coordmath import get_locus_string
from prymer.api.coordmath import require_same_refname


def test_get_locus_string() -> None:
    intervals = [Interval("chr1", 1, 2), Interval("chr2", 3, 4)]
    assert [get_locus_string(i) for i in intervals] == ["chr1:1-2", "chr2:3-4"]


def test_assert_same_refname_works() -> None:
    intervals = [Interval("chr1", 1, 2), Interval("chr1", 3, 4)]
    require_same_refname(*intervals)


def test_assert_same_refname_works_neg() -> None:
    intervals = [Interval("chr1", 1, 2), Interval("chr2", 3, 4)]
    with pytest.raises(ValueError):
        require_same_refname(*intervals)
