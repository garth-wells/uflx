"""Test finite element."""

import pytest
from utils import point, interval, triangle, tetrahedron, quadrilateral, LagrangeElement, hexahedron


@pytest.mark.parametrize(
    "cell",
    [
        point,
        interval,
        triangle,
        quadrilateral,
        tetrahedron,
        hexahedron,
    ],
)
def test_lagrange_element(cell):
    element = LagrangeElement(cell, 2)

    assert element.cell == cell
    assert element.reference_value_shape == ()
    assert element.embedded_lagrange_superdegree == 2
