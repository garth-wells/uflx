"""Test solving a full Poisson problem."""

import numpy as np
import pytest
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

from basix_uflx import element


def apply_bcs(matrix, vector, bcs):
    """Apply boundary conditions."""
    ndofs = matrix.shape[0]
    for i in range(ndofs):
        if i in bcs:
            matrix[i, :] = 0.0
            matrix[i, i] = 1.0
            vector[i] = bcs[i]
        else:
            for dof, value in bcs.items():
                vector[i] -= value * matrix[i, dof]
                matrix[i, dof] = 0.0
    return matrix, vector


@pytest.mark.parametrize("npoints", [1, 3, 8])
@pytest.mark.parametrize("degree", range(1, 5))
def test_poisson_problem_square(
    npoints, degree, code_dir, create_mesh, create_function_space, assemble_matrix, assemble_vector
):
    """Test a full Poisson problem solve.

    This tests a manufactured problem with the solution u = (x - y) ** degree.
    For this problem, Δu = 2 * degree * (degree - 1) * (x - y) ** (degree - 2) and
    we use Dirichlet BCs on all four sides of a unit square.
    """
    points = np.array(
        [[i / npoints, j / npoints] for j in range(npoints + 1) for i in range(npoints + 1)]
    )

    cells = []
    for j in range(npoints):
        for i in range(npoints):
            origin = j * (npoints + 1) + i
            cells.append([origin, origin + 1, origin + npoints + 2])
            cells.append([origin, origin + npoints + 2, origin + npoints + 1])

    e = element("Lagrange", "triangle", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "triangle", 1, shape=(2,)))

    mesh = create_mesh(
        points=points,
        cells=cells,
        domain=domain,
    )

    dofmap = {"point": {(i,): [i] for i, _ in enumerate(points)}}
    dof_locations = {i: p for i, p in enumerate(points)}
    dof_n = points.shape[0]
    if degree > 1:
        ndofs_interval = degree - 1
        dofmap["interval"] = {}
        for cell in cells:
            for edge in [(cell[0], cell[1]), (cell[0], cell[2]), (cell[1], cell[2])]:
                if edge not in dofmap["interval"] and edge[::-1] not in dofmap["interval"]:
                    dofmap["interval"][edge] = [dof_n + i for i in range(ndofs_interval)]
                    a = points[edge[0]]
                    b = points[edge[1]]
                    for i in range(ndofs_interval):
                        dof_locations[dof_n + i] = a + (i + 1) * (b - a) / (ndofs_interval + 1)
                    dof_n += ndofs_interval
    if degree > 2:
        ndofs_triangle = (degree - 2) * (degree - 1) // 2
        dofmap["triangle"] = {}
        for cell in cells:
            dofmap["triangle"][tuple(cell)] = [dof_n + i for i in range(ndofs_triangle)]
            a = points[cell[0]]
            b = points[cell[1]]
            c = points[cell[2]]
            for j in range(1, degree - 1):
                for i in range(1, degree - j):
                    dof_locations[dof_n] = a + i * (b - a) / degree + j * (c - a) / degree
                    dof_n += 1

    ndofs = (npoints * degree + 1) ** 2
    assert dof_n == ndofs

    space = function_space(domain, e)

    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(grad(u), grad(v)) * dx
    x = SpatialCoordinate(2)
    rhs = -2 * degree * (degree - 1) * (x[0] - x[1]) ** (degree - 2) * v * dx

    wrapped_space = create_function_space(space=space, mesh=mesh, dofmap=dofmap)

    matrix = assemble_matrix(form, wrapped_space, wrapped_space, code_dir)
    vector = assemble_vector(rhs, wrapped_space, code_dir)

    matrix, vector = apply_bcs(
        matrix,
        vector,
        {
            d: (p[0] - p[1]) ** degree
            for d, p in dof_locations.items()
            if np.isclose(p[0], 0.0)
            or np.isclose(p[0], 1.0)
            or np.isclose(p[1], 0.0)
            or np.isclose(p[1], 1.0)
        },
    )

    solution = np.linalg.inv(matrix) @ vector

    for d, (x, y) in dof_locations.items():
        assert np.isclose(solution[d], (x - y) ** degree)
