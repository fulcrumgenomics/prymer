from prymer.api.coordmath import get_closed_end


def test_get_closed_end() -> None:
    assert get_closed_end(start=10, length=0) == 9
    assert get_closed_end(start=10, length=1) == 10
    assert get_closed_end(start=10, length=10) == 19
