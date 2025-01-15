from typing import TypeVar

import pytest

from prymer import MinOptMax

Numeric = TypeVar("Numeric", int, float)


@pytest.mark.parametrize(
    "valid_min,valid_opt,valid_max,expected_obj",
    [
        (0.0, 1.0, 2.0, MinOptMax(0.0, 1.0, 2.0)),  # min < opt < max
        (0.0, 0.0, 2.0, MinOptMax(0.0, 0.0, 2.0)),  # min == opt < max
        (1.0, 1.0, 2.0, MinOptMax(1.0, 1.0, 2.0)),  # min == opt < max non-zero
        (1.0, 1.0, 1.0, MinOptMax(1.0, 1.0, 1.0)),  # min == opt == max
        (-100, -5, 42, MinOptMax(-100, -5, 42)),  # min < opt < max (negative ints)
    ],
)
def test_minmaxopt_valid(
    valid_min: Numeric, valid_opt: Numeric, valid_max: Numeric, expected_obj: MinOptMax[Numeric]
) -> None:
    """Test MinOptMax construction with valid input"""
    test_obj = MinOptMax(valid_min, valid_opt, valid_max)
    assert test_obj == expected_obj


@pytest.mark.parametrize(
    "invalid_min,invalid_opt, invalid_max",
    [
        (3.0, 12.0, 10.0),  # opt > max
        (10.0, 5.0, 1.0),  # min > max
        (1, 50.0, "string_value"),  # mixed types
        (1.0, 0.0, 2.0),  # opt is 0.0 and 0.0<min<max
    ],
)
def test_minoptmax_raises(invalid_min: Numeric, invalid_opt: Numeric, invalid_max: Numeric) -> None:
    """Test that MinMaxOpt constructor raises an error with invalid input"""
    with pytest.raises((ValueError, TypeError)):
        MinOptMax(min=invalid_min, opt=invalid_opt, max=invalid_max)


def test_str_repr() -> None:
    test_obj = MinOptMax(min=1, opt=2, max=3)
    assert test_obj.__str__() == "(min:1, opt:2, max:3)"


def test_iter_repr() -> None:
    test_obj = MinOptMax(min=1, opt=2, max=3)
    test_iter = test_obj.__iter__()
    assert next(test_iter) == 1
    assert next(test_iter) == 2
    assert next(test_iter) == 3
    with pytest.raises(StopIteration):
        next(test_iter)
