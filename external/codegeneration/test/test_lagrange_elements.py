"""Test implementation of Lagrange elements used in tests."""

import numpy as np
import pytest
from cffi import FFI
from uflx import (
    SpatialCoordinate,
    TestFunction,
    TrialFunction,
    coordinate_element,
    dx,
    function_space,
    grad,
    inner,
)

from uflx_codegeneration.utils import index


@pytest.mark.parametrize("degree", range(1, 4))
def test_dofs(lagrange_element, degree):
    """Test that basis functions are 0 and 1 at DOF points."""
    element = lagrange_element("triangle", degree)

    points = np.zeros((element.dim, 2))

    if degree > 0:
        points[0] = [0.0, 0.0]
        points[1] = [1.0, 0.0]
        points[2] = [0.0, 1.0]

    n = 3
    if degree > 1:
        for i in range(1, degree):
            points[n] = [i / degree, 0.0]
            n += 1
        for i in range(1, degree):
            points[n] = [0.0, i / degree]
            n += 1
        for i in range(1, degree):
            points[n] = [1.0 - i / degree, i / degree]
            n += 1
    if degree > 2:
        for j in range(1, degree - 1):
            for i in range(1, degree - j):
                points[n] = [i / degree, j / degree]
                n += 1

    assert n == element.dim

    table = element.tabulate(0, points)
    assert np.allclose(table[0, :, :, 0], np.eye(element.dim, element.dim))


@pytest.mark.parametrize("nderivs", range(1, 6))
@pytest.mark.parametrize("degree", range(1, 4))
def test_derivatives(lagrange_element, degree, nderivs):
    """Test that basis functions are 0 and 1 at DOF points."""
    element = lagrange_element("triangle", degree)

    n = 15
    points = np.array([[i / n, j / n] for j in range(n + 1) for i in range(n + 1 - j)])
    eps = 1e-8
    points_x = np.array([[i + eps, j] for i, j in points])
    points_y = np.array([[i, j + eps] for i, j in points])

    table = element.tabulate(nderivs, points)
    table_x = element.tabulate(nderivs, points_x)
    table_y = element.tabulate(nderivs, points_y)

    # test d/dx
    for ix in range(1, nderivs + 1):
        iy = nderivs - ix
        assert np.allclose(
            table[index(ix, iy)],
            (table_x[index(ix - 1, iy)] - table[index(ix - 1, iy)]) / eps,
            atol=eps * 100,
        )

    # test d/dy
    for iy in range(1, nderivs + 1):
        ix = nderivs - iy
        assert np.allclose(
            table[index(ix, iy)],
            (table_y[index(ix, iy - 1)] - table[index(ix, iy - 1)]) / eps,
            atol=eps * 100,
        )
