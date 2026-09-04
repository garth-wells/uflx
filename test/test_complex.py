"""Test complex numbers."""
import pytest

from uflx.expressions import Integer, ComplexScalar
from uflx.complex import take_real_part, take_imaginary_part
from uflx.graphs import generate_graph

@pytest.fixture
def z():
    a = Integer(4)
    b = Integer(6)
    return generate_graph(ComplexScalar(a, b))


def test_real_part(z):
    """Test taking real part."""
    assert take_real_part(z).root == 4


def test_imaginary_part(z):
    """Test taking imaginary part."""
    assert take_imaginary_part(z).root == 6
