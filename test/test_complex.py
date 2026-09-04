"""Test complex numbers."""

import pytest

from uflx.complex import take_imaginary_part, take_real_part
from uflx.expressions import ComplexScalar, Integer
from uflx.graphs import generate_graph


@pytest.fixture
def z():
    """The complex number 4+6j."""
    a = Integer(4)
    b = Integer(6)
    return generate_graph(ComplexScalar(a, b))


def test_real_part(z):
    """Test taking real part."""
    assert take_real_part(z).root == 4


def test_imaginary_part(z):
    """Test taking imaginary part."""
    assert take_imaginary_part(z).root == 6
