"""Test complex values."""

import numpy as np

from uflx.complex import take_imaginary_part, take_real_part
from uflx.expressions import ComplexScalar, Integer, to_scalar


def test_real_part():
    """Test taking real part."""
    assert take_real_part(ComplexScalar(Integer(4), Integer(6))) == 4


def test_imaginary_part():
    """Test taking imaginary part."""
    assert take_imaginary_part(ComplexScalar(Integer(4), Integer(6))) == 6


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
