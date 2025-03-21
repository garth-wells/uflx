"""Test function space."""

import pytest
from utils import triangle, quadrilateral, tetrahedron, LagrangeElement
from uflx import FunctionSpace, EmbeddedCell

@pytest.mark.parametrize(("cell", "gdim"), [
    (triangle, 2),
    (triangle, 3),
    (quadrilateral, 2),
    (quadrilateral, 3),
    (tetrahedron, 3),
])
def test_function_space(cell, gdim):
    space = FunctionSpace([(EmbeddedCell(cell, gdim), LagrangeElement(cell, 1))])
    assert len(space.elements) == len(space.domains) == 1

@pytest.mark.parametrize("gdim", [2, 3])
def test_function_space_multiple_cells(gdim):
    space = FunctionSpace([
        (EmbeddedCell(triangle, gdim), LagrangeElement(triangle, 2)),
        (EmbeddedCell(quadrilateral, gdim), LagrangeElement(quadrilateral, 2)),
    ])
    assert len(space.elements) == len(space.domains) == 2

