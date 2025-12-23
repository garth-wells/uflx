"""Test cell."""

import pytest
from utils import interval, point, tetrahedron, triangle


@pytest.mark.parametrize(
    "cell",
    [
        point,
        interval,
        triangle,
        tetrahedron,
    ],
)
def test_euler_characteristic(cell):
    """Test Euler characteristic of cell."""
    match cell.topological_dimension:
        case 0:
            assert len(cell.vertices) == 1
        case 1:
            assert len(cell.vertices) - len(cell.edges) == 1
        case 2:
            assert len(cell.vertices) - len(cell.edges) + len(cell.faces) == 1
        case 3:
            assert len(cell.vertices) - len(cell.edges) + len(cell.faces) == 2
