import numpy as np
import uflx
from uflx.expressions import to_scalar


def test_complex_scalar_add():
    five = to_scalar(5)

    z = five + 2j
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), 2)

    z = 2j + five
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), 2)


def test_complex_scalar_sub():
    five = to_scalar(5)

    z = five - 2j
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), -2)

    z = 2j - five
    assert np.isclose(z.re.as_float(), -5)
    assert np.isclose(z.im.as_float(), 2)


def test_complex_scalar_mult():
    five = to_scalar(5)

    z = five * 2j
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 10)

    z = 2j * five
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 10)


def test_complex_scalar_div():
    five = to_scalar(5)

    z = five / 2j
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), -2.5)

    z = 2j / five
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 0.4)


def test_complex_scalar_neg():
    z = to_scalar(3 - 2j)

    assert np.isclose(z.re.as_float(), 3)
    assert np.isclose(z.im.as_float(), -2)

    assert np.isclose((-z).re.as_float(), -3)
    assert np.isclose((-z).im.as_float(), 2)
