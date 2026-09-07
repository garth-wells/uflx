"""Test geometry."""

import pytest

from uflx import coordinate_element
from uflx.expressions import RealScalar
from uflx.geometry import Jacobian, JacobianInverse, JacobianInverseTranspose, JacobianTranspose
from uflx.points import Point

cells_and_gdims = [
    ("interval", 1),
    ("interval", 2),
    ("interval", 3),
    ("interval", 4),
    ("triangle", 2),
    ("triangle", 3),
    ("triangle", 4),
    ("quadrilateral", 2),
    ("quadrilateral", 3),
    ("quadrilateral", 4),
    ("tetrahedron", 3),
    ("tetrahedron", 4),
    ("hexahedron", 3),
    ("hexahedron", 4),
]


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * gdim)
    j = Jacobian(domain, point)
    assert j.value_shape == (gdim, tdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (gdim, tdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_inverse_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * gdim)
    j = JacobianInverse(domain, point)
    assert j.value_shape == (tdim, gdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (tdim, gdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_tranpose_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse transpose."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * gdim)
    j = JacobianTranspose(domain, point)
    assert j.value_shape == (tdim, gdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (tdim, gdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_inverse_transpose_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse transpose."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * gdim)
    j = JacobianInverseTranspose(domain, point)
    assert j.value_shape == (gdim, tdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (gdim, tdim)
