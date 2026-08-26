import os
import random
import pytest
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
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
from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractReferenceMappedFiniteElement
from uflx.function_spaces import AbstractReferenceMappedFunctionSpace
from uflx.integrals import AbstractIntegral
from uflx.domains import AbstractCoordinateElement
import uflx_codegeneration
import hashlib

from basix_uflx import element


class Mesh(NamedTuple):
    """A mesh."""
    points: npt.NDArray[float]
    cells: list[list[int]]
    domain: AbstractCoordinateElement


class FunctionSpace(NamedTuple):
    """A function space."""
    mesh: Mesh
    space: AbstractReferenceMappedFunctionSpace
    dofmap: dict[str, dict[tuple[int, ...], list[int]]]


def assemble_code(form, code_dir, filename=None):
    ffi = FFI()
    code, signatures = uflx_codegeneration.generate(form)

    if filename is None:
        h = hashlib.sha1(code.encode("utf-8"))
        filename = f"basix_uflx_test_{h.hexdigest()}"
        while os.path.isfile(os.path.join(code_dir, f"{filename}.c")):
            filename += "_"

    ffi.cdef("\n".join(signatures.values()))
    ffi.set_source(filename, code)
    so = ffi.compile(code_dir)

    lib = ffi.dlopen(so)
    empty = np.zeros(0)

    def tabulate(local_tensor, coords):
        local_tensor[:] = 0.0
        lib.tabulate_tensor_f64(
            ffi.cast("double*", local_tensor.ctypes.data),
            ffi.cast("double*", empty.ctypes.data),
            ffi.cast("double*", empty.ctypes.data),
            ffi.cast("double*", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )

    return tabulate


def assemble_matrix(
    form: AbstractIntegral,
    test_function_space: FunctionSpace,
    trial_function_space: FunctionSpace,
    code_dir: str,
) -> npt.NDArray[float]:
    tabulate = assemble_code(form, code_dir)

    mesh = test_function_space.mesh
    assert trial_function_space.mesh == mesh

    test_dim = test_function_space.space.elements[0].dim
    trial_dim = trial_function_space.space.elements[0].dim
    assert test_function_space.space.domain == trial_function_space.space.domain == mesh.domain
    gdim = mesh.domain.geometric_dimension
    assert len(mesh.domain.cells) == 1
    mesh_cell = mesh.domain.cells[0]
    tdim = mesh_cell.topological_dimension

    test_ndofs = sum(len(entity_dofs) for dofs in test_function_space.dofmap.values() for entity_dofs in dofs.values())
    trial_ndofs = sum(len(entity_dofs) for dofs in trial_function_space.dofmap.values() for entity_dofs in dofs.values())

    npoints = mesh.domain.cells[0].sub_entity_count(0)

    local_mat = np.zeros((test_dim, trial_dim))
    coords = np.zeros((npoints, tdim))

    matrix = np.zeros((test_ndofs, trial_ndofs))

    for cell in mesh.cells:
        for i, j in enumerate(cell):
            for k, p in enumerate(mesh.points[j]):
                coords[i, k] = p
        test_dofs = []
        trial_dofs = []
        for d in range(tdim + 1):
            for et, vs in zip(mesh_cell.sub_entities(d), mesh_cell.sub_entity_vertices(d)):
                cell_vs = tuple(cell[i] for i in vs)
                test_dofs += test_function_space.dofmap.get(et.name, {}).get(cell_vs, [])
                trial_dofs += trial_function_space.dofmap.get(et.name, {}).get(cell_vs, [])

        tabulate(local_mat, coords)
        for test_dof, row in zip(test_dofs, local_mat):
            for trial_dof, entry in zip(trial_dofs, row):
                matrix[test_dof, trial_dof] += entry

    return matrix


def assemble_vector(
    form: AbstractIntegral,
    test_function_space: FunctionSpace,
    code_dir: str,
) -> npt.NDArray[float]:
    tabulate = assemble_code(form, code_dir)

    mesh = test_function_space.mesh

    test_dim = test_function_space.space.elements[0].dim
    assert test_function_space.space.domain == mesh.domain
    gdim = mesh.domain.geometric_dimension
    assert len(mesh.domain.cells) == 1
    mesh_cell = mesh.domain.cells[0]
    tdim = mesh_cell.topological_dimension

    test_ndofs = sum(len(entity_dofs) for dofs in test_function_space.dofmap.values() for entity_dofs in dofs.values())

    npoints = mesh.domain.cells[0].sub_entity_count(0)

    local_vec = np.zeros(test_dim)
    coords = np.zeros((npoints, tdim))

    vector = np.zeros(test_ndofs)

    for cell in mesh.cells:
        for i, j in enumerate(cell):
            for k, p in enumerate(mesh.points[j]):
                coords[i, k] = p
        test_dofs = []
        for d in range(tdim + 1):
            for et, vs in zip(mesh_cell.sub_entities(d), mesh_cell.sub_entity_vertices(d)):
                cell_vs = tuple(cell[i] for i in vs)
                test_dofs += test_function_space.dofmap.get(et.name, {}).get(cell_vs, [])

        tabulate(local_vec, coords)
        for test_dof, entry in zip(test_dofs, local_vec):
            vector[test_dof] += entry

    return vector


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


@pytest.mark.parametrize("npoints", [1, 3, 8])
@pytest.mark.parametrize("degree", range(1, 3))
def test_poisson_problem_square(npoints, degree, code_dir):
    """Test a full Poisson problem solve.

    This tests a manufactured problem with the solution u = (x - y) ** degree.
    For this problem, Δu = 2 * degree * (degree - 1) * (x - y) ** (degree - 2) and
    we use Dirichlet BCs on all four sides of a unit square.
    """

    points = np.array([[i / npoints, j / npoints] for j in range(npoints+1) for i in range(npoints+1)])

    cells = []
    for j in range(npoints):
        for i in range(npoints):
            origin = j * (npoints + 1) + i
            cells.append([origin, origin + 1, origin + npoints + 2])
            cells.append([origin, origin + npoints + 2, origin + npoints + 1])
    boundary_points = [
        j * (npoints + 1) + i
        for j in range(npoints + 1)
        for i in (range(npoints+1) if j in [0, npoints] else [0, npoints])
    ]

    e = element("Lagrange", "triangle", degree)
    domain = coordinate_element(element("Lagrange", "triangle", 1, shape=(2,)))

    mesh = Mesh(
        points=points,
        cells=cells,
        domain=domain,
    )

    dofmap = {"point": {(i, ): [i] for i in range(ndofs)}}
    dof_n = (npoints + 1) ** 2
    if degree > 1:
        dofmap["interval"] = {}
        for cell in cells:
            for edge in [(cell[0], cell[1]), (cell[0], cell[2]), (cell[1], cell[2])]:
                if edge not in dofmap["interval"] and edge[::-1] not in dofmap["interval"]:
                    dof_n +++++

    ndofs = (npoints + 1) ** 2


    space = function_space(domain, e)

    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(grad(u), grad(v)) * dx
    x = SpatialCoordinate(2)
    rhs = 2 * degree * (degree - 1) * (x[0] - x[1]) ** (degree - 2) * v * dx

    wrapped_space = FunctionSpace(space=space, mesh=mesh, dofmap=dofmap)

    matrix = assemble_matrix(form, wrapped_space, wrapped_space, code_dir)
    vector = assemble_vector(rhs, wrapped_space, code_dir)

    matrix, vector = apply_bcs(matrix, vector, {
        d: (points[i][0] - points[i][1]) ** degree
        for i in boundary_points
        for d in dofmap["point"][(i,)]
    })

    solution = np.linalg.inv(matrix) @ vector

    for i in range(1, npoints):
        for j in range(1, npoints):
            index = (npoints + 1) * j + i
            x, y = points[index]
            for d in dofmap["point"][(index,)]:
                assert np.isclose(solution[d], (x - y) ** degree)
