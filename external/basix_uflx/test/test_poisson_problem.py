import random
import pytest

import numpy as np
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
import uflx_codegeneration

from basix_uflx import element


def assemble_code(form, filename, code_dir):
    ffi = FFI()
    code, signatures = uflx_codegeneration.generate(form)

    ffi.cdef("\n".join(signatures.values()))
    ffi.set_source(filename, code)
    so = ffi.compile(code_dir)

    lib = ffi.dlopen(so)
    empty = np.zeros(0)

    def tabulate(local_mat, coords):
        local_mat[:] = 0.0
        lib.tabulate_tensor_f64(
            ffi.cast("double*", local_mat.ctypes.data),
            ffi.cast("double*", empty.ctypes.data),
            ffi.cast("double*", empty.ctypes.data),
            ffi.cast("double*", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )

    return tabulate

def apply_bcs(matrix, vector, bcs):
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


@pytest.mark.parametrize("n", [1, 3, 8, 15])
def test_poisson_problem_square(n, code_dir):
    """Test a full Poisson problem solve.

    This tests a manufactured problem with the solution u = x - y.
    For this problem, Δu = 0 and we use Dirichlet BCs on all four sides of a unit square.
    """

    points = np.array([[i / n, j / n] for j in range(n+1) for i in range(n+1)])

    cells = []
    for j in range(n):
        for i in range(n):
            origin = j * (n + 1) + i
            cells.append([origin, origin + 1, origin + n + 2])
            cells.append([origin, origin + n + 2, origin + n + 1])
    boundary_points = [
        j * (n + 1) + i
        for j in range(n + 1)
        for i in (range(n+1) if j in [0, n] else [0, n])
    ]

    e = element("Lagrange", "triangle", 1)
    ndofs = (n + 1) ** 2
    dofmap = {"vertices": [[i] for i in range(ndofs)]}
    # random.shuffle(dofmap["vertices"])  # Use randomly numbered DOFs

    domain = coordinate_element(element("Lagrange", "triangle", 1, shape=(2,)))
    space = function_space(domain, e)

    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(grad(u), grad(v)) * dx
    rhs = 0 * v * dx

    mat_tabulate = assemble_code(form, f"test_poisson_problem_square_{n}_mat", code_dir)
    vec_tabulate = assemble_code(rhs, f"test_poisson_problem_square_{n}_vec", code_dir)


    local_mat = np.zeros((3, 3))
    local_vec = np.zeros(3)
    coords = np.zeros((3, 2))

    matrix = np.zeros((ndofs, ndofs))
    vector = np.zeros(ndofs)

    for cell in cells:
        for i, j in enumerate(cell):
            for k, p in enumerate(points[j]):
                coords[i, k] = p
        dofs = []
        for v in cell:
            dofs += dofmap["vertices"][v]

        mat_tabulate(local_mat, coords)
        print(local_mat)
        for row_dof, row in zip(dofs, local_mat):
            for col_dof, entry in zip(dofs, row):
                matrix[row_dof, col_dof] += entry

        vec_tabulate(local_vec, coords)
        for row_dof, entry in zip(dofs, local_vec):
            vector[row_dof] += entry

    print(matrix)
    print(vector)

    matrix, vector = apply_bcs(matrix, vector, {
        d: points[i][0] - points[i][1]
        for i in boundary_points
        for d in dofmap["vertices"][i]
    })

    solution = np.linalg.inv(matrix) @ vector

    for i in range(1, n):
        for j in range(1, n):
            index = (n + 1) * j + i
            x, y = points[index]
            for d in dofmap["vertices"][index]:
                print((float(x), float(y)), solution[d], x - y)
    print()

    for i in range(1, n):
        for j in range(1, n):
            index = (n + 1) * j + i
            x, y = points[index]
            for d in dofmap["vertices"][index]:
                assert np.isclose(solution[d], x - y)
