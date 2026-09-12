"""Check higher-order demo DOF continuity independently of kernel assembly."""

import importlib
from pathlib import Path

import basix
import numpy as np
import pytest


@pytest.fixture
def builders(monkeypatch):
    """Import the loose demo modules without permanently changing the search path."""
    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "demo"))
    mesh = importlib.import_module("assemble_mesh_gpu")
    nodal = importlib.import_module("lagrange_dofmap")
    return mesh, nodal


@pytest.mark.parametrize("degree", [1, 2, 3, 4])
@pytest.mark.parametrize("permute", [False, True])
def test_nodal_positions_and_shared_ids(builders, degree, permute):
    """Physical node positions identify exactly the same shared DOFs as topology."""
    mesh, nodal = builders
    points, cells = mesh.build_mesh(2)
    cells = np.asarray(cells)
    if permute:
        rng = np.random.default_rng(513)
        cells = np.stack([rng.permutation(cell) for cell in cells])
    flat, ndofs, nglobal = nodal.build_dofmap(cells.tolist(), len(points), degree)
    assert nglobal == (2 * degree + 1) ** 3
    assert ndofs == (degree + 1) * (degree + 2) * (degree + 3) // 6
    assert flat.dtype == np.int32
    element = basix.create_element(
        basix.ElementFamily.P,
        basix.CellType.tetrahedron,
        degree,
        basix.LagrangeVariant.equispaced,
    )
    # Verify Basix's basis ordering at its interpolation points, then independently
    # interpolate those points into physical space for every cell.
    np.testing.assert_allclose(
        element.tabulate(0, element.points)[0, :, :, 0], np.eye(ndofs), atol=1e-12
    )
    barycentric = np.column_stack((1 - element.points.sum(axis=1), element.points))
    xyz = np.einsum("iv,cvd->cid", barycentric, points[cells])
    lattice = np.rint(xyz * (2 * degree)).astype(np.int64)
    np.testing.assert_allclose(xyz * (2 * degree), lattice, atol=1e-12)
    grid_size = 2 * degree + 1
    keys = ((lattice[..., 0] * grid_size + lattice[..., 1]) * grid_size + lattice[..., 2]).reshape(
        -1
    )
    pairs = np.unique(np.column_stack((flat, keys)), axis=0)
    # Each topological ID maps to exactly one physical node, and vice versa.
    assert len(pairs) == nglobal
    assert len(np.unique(pairs[:, 0])) == nglobal
    assert len(np.unique(pairs[:, 1])) == nglobal


@pytest.mark.parametrize("degree", [1, 2])
def test_preserves_p1_p2_numbering(builders, degree):
    """Retain the existing low-order map and global vertex IDs."""
    mesh, nodal = builders
    points, cells = mesh.build_mesh(2)
    old, old_dim, old_size = mesh.build_dofmap(cells, len(points), degree)
    new, new_dim, new_size = nodal.build_dofmap(cells, len(points), degree)
    np.testing.assert_array_equal(new, old)
    assert (new_dim, new_size) == (old_dim, old_size)
