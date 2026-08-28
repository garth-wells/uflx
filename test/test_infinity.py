import pytest
from uflx.utils import infinity, Infinity, negative_infinity, NegativeInfinity


def test_is():
    assert infinity is infinity
    assert infinity is Infinity()
    assert negative_infinity is negative_infinity
    assert negative_infinity is NegativeInfinity()

    assert negative_infinity is not infinity
    assert negative_infinity is not Infinity()
    assert infinity is not negative_infinity
    assert infinity is not NegativeInfinity()


def test_add():
    assert infinity + -1 is infinity
    assert infinity + 0 is infinity
    assert infinity + 1 is infinity
    assert infinity + 2 is infinity
    assert infinity + 3 is infinity
    assert infinity + infinity is infinity

    with pytest.raises(TypeError):
        infinity + None

    assert negative_infinity + -1 is negative_infinity
    assert negative_infinity + 0 is negative_infinity
    assert negative_infinity + 1 is negative_infinity
    assert negative_infinity + 2 is negative_infinity
    assert negative_infinity + 3 is negative_infinity
    assert negative_infinity + negative_infinity is negative_infinity

    with pytest.raises(TypeError):
        negative_infinity + None

    with pytest.raises(ValueError):
        negative_infinity + infinity


def test_subtract():
    assert infinity - -1 is infinity
    assert infinity - 0 is infinity
    assert infinity - 1 is infinity
    assert infinity - 2 is infinity
    assert infinity - 3 is infinity
    assert infinity - negative_infinity is infinity

    with pytest.raises(TypeError):
        infinity - None

    assert negative_infinity - -1 is negative_infinity
    assert negative_infinity - 0 is negative_infinity
    assert negative_infinity - 1 is negative_infinity
    assert negative_infinity - 2 is negative_infinity
    assert negative_infinity - 3 is negative_infinity
    assert negative_infinity - infinity is negative_infinity

    with pytest.raises(TypeError):
        negative_infinity - None

    with pytest.raises(ValueError):
        negative_infinity - negative_infinity
    with pytest.raises(ValueError):
        infinity - infinity
