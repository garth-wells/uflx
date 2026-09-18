"""Conforming P1-P4 equispaced tetrahedral DOF maps for the assembly demos."""

import basix
import numpy as np


def build_dofmap(
    cells: list[tuple[int, int, int, int]], nverts: int, degree: int
) -> tuple[np.ndarray, int, int]:
    """Number Basix nodes by their barycentric weights on global mesh vertices.

    Integer barycentric weights identify edge, face and cell-interior nodes
    exactly for equispaced elements. Sort vertex/weight pairs together to
    identify the same node across differently oriented neighbouring cells.
    This handles reversed edge nodes and permuted face nodes without assuming
    a particular local edge or face ordering. Vertex DOFs retain vertex IDs;
    remaining DOFs are numbered by entity dimension.

    Args:
        cells: Tetrahedral connectivity using global vertex IDs.
        nverts: Number of mesh vertices.
        degree: Lagrange degree, from 1 to 4.

    Returns:
        Flattened int32 cell-to-DOF map, local dimension and global DOF count.
    """
    if degree not in (1, 2, 3, 4):
        raise ValueError("degree must be 1, 2, 3 or 4")
    element = basix.create_element(
        basix.ElementFamily.P,
        basix.CellType.tetrahedron,
        degree,
        basix.LagrangeVariant.equispaced,
    )
    # For scalar nodal Lagrange elements the interpolation points follow the
    # basis-function/DOF ordering. Check this contract rather than assuming it.
    if not element.interpolation_is_identity or element.points.shape != (element.dim, 3):
        raise ValueError("Expected a scalar nodal element with identity interpolation")
    points = element.points
    barycentric = np.column_stack((1 - points.sum(axis=1), points))
    weights = np.rint(degree * barycentric).astype(np.int64)
    if not np.allclose(weights, degree * barycentric, rtol=0, atol=1e-12):
        raise ValueError("Interpolation points are not on the equispaced barycentric lattice")
    connectivity = np.asarray(cells, dtype=np.int64).reshape(-1, 4)
    if nverts < 0 or (
        connectivity.size and (connectivity.min() < 0 or connectivity.max() >= nverts)
    ):
        raise ValueError("Invalid vertex count or connectivity")
    if connectivity.size and np.any(np.diff(np.sort(connectivity, axis=1), axis=1) == 0):
        raise ValueError("A tetrahedron must have four distinct vertices")
    ncells = len(connectivity)
    dofs = np.empty((ncells, element.dim), dtype=np.int32)
    support_sizes = np.count_nonzero(weights, axis=1)
    for local in np.flatnonzero(support_sizes == 1):
        vertex = int(np.flatnonzero(weights[local])[0])
        dofs[:, local] = connectivity[:, vertex]
    nglobal = nverts
    for support_size in (2, 3, 4):
        local_dofs = np.flatnonzero(support_sizes == support_size)
        if not len(local_dofs) or not ncells:
            continue
        keys = np.empty((ncells, len(local_dofs), 2 * support_size), dtype=np.int64)
        for slot, local in enumerate(local_dofs):
            support = np.flatnonzero(weights[local])
            vertices = connectivity[:, support]
            order = np.argsort(vertices, axis=1)
            keys[:, slot, :support_size] = np.take_along_axis(vertices, order, axis=1)
            keys[:, slot, support_size:] = weights[local, support][order]
        unique, inverse = np.unique(keys.reshape(-1, 2 * support_size), axis=0, return_inverse=True)
        if nglobal + len(unique) > np.iinfo(np.int32).max:
            raise ValueError("DOF map exceeds int32 indexing capacity")
        dofs[:, local_dofs] = nglobal + inverse.reshape(ncells, len(local_dofs))
        nglobal += len(unique)
    return dofs.reshape(-1), element.dim, nglobal
