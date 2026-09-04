"""Test points."""

import pytest

from uflx.expressions import Integer
from uflx.points import RD, Point


@pytest.mark.parametrize("dim", range(5))
def test_rd(dim):
    """Test R^d set of points."""
    points = RD(dim)
    assert points.geometric_dimension == dim
    assert points.npoints == "Infinity"


@pytest.mark.parametrize("dim", range(5))
def test_point(dim):
    """Test a point."""
    point = Point([Integer(i) for i in range(dim)])

    assert point.dim == dim
    assert point.points_set.geometric_dimension == dim
    assert point.points_set.npoints == "Infinity"
