"""Test complex values."""

import numpy as np

from uflx.expressions import to_scalar


def test_complex_scalar_add():
    """Test addition with complex scalars."""
    five = to_scalar(5)

    z = five + 2j
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), 2)

    z = 2j + five
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), 2)


def test_complex_scalar_sub():
    """Test subtraction with complex scalars."""
    five = to_scalar(5)

    z = five - 2j
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), -2)

    z = 2j - five
    assert np.isclose(z.re.as_float(), -5)
    assert np.isclose(z.im.as_float(), 2)


def test_complex_scalar_mult():
    """Test multiplication with complex scalars."""
    five = to_scalar(5)

    z = five * 2j
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 10)

    z = 2j * five
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 10)


def test_complex_scalar_div():
    """Test division with complex scalars."""
    five = to_scalar(5)

    z = five / 2j
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), -2.5)

    z = 2j / five
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 0.4)


def test_complex_scalar_neg():
    """Test negation with complex scalars."""
    z = to_scalar(3 - 2j)

    assert np.isclose(z.re.as_float(), 3)
    assert np.isclose(z.im.as_float(), -2)

    assert np.isclose((-z).re.as_float(), -3)
    assert np.isclose((-z).im.as_float(), 2)
