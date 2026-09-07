"""Assemble a real global CSR stiffness matrix using uflx_mlir.gpu_assembly's
single-call, whole-mesh CSR assembly kernel (generate_csr_assembly_module),
on a structured tetrahedral box mesh at an arbitrary supported Lagrange
degree (P2 by default).

Context: gpu_assembly.py's own tests (test_gpu_assembly.py) exercise
generate_csr_assembly_module (the kernel used here) on only a genuine but
tiny two-cell mesh (see
test_generate_csr_assembly_module_matches_quadrature_reference_two_cells)
-- enough to prove that dofs shared between cells (each contributing an
atomic add to the same matrix entry) combine correctly, but not at any
real mesh scale. This script is that next step: a real (if modest-sized)
mesh, a real P1/P2 dof map with genuinely shared dofs, and a real global
CSR matrix built with a SINGLE call to generate_csr_assembly_module's
kernel -- which itself loops over every cell, refreshes that cell's own
geometry, accumulates its local stiffness block, and flushes it into the
CSR matrix via the same binary-search-and-atomic-add scatter
generate_csr_entry_module's single-cell kernel uses (see that function's
docstring). There is no host-side Python loop over cells or dof pairs at
all -- see "Why not the real GPU path" below for what's still missing
compared to running this on an actual GPU.

Mesh: a unit cube split into n x n x n sub-cubes, each split into 6
tetrahedra via the standard Kuhn/Freudenthal triangulation (all 6 sharing
the cube's main space diagonal). Orientation is not controlled for --
irrelevant here, since uflx_mlir.geometry's affine-Poisson geometry
extraction always takes abs(detJ) (see demo/README.md's "What this does
and doesn't cover yet").

Dof numbering: P1 uses vertex dofs only (global dof id == vertex id). P2
adds one dof per mesh edge (shared between however many cells touch that
edge), numbered after the vertex dofs; a fresh global id is assigned the
first time each edge is seen, keyed by its sorted pair of global vertex
ids so the SAME id is reused by every cell touching that edge. The
per-cell local-dof-slot -> local-edge mapping below (LOCAL_EDGES) was
cross-checked against basix directly (both
basix.cell.sub_entity_connectivity(tetrahedron)[1] and
basix.create_element(..., degree=2, equispaced).entity_dofs), not
guessed -- local dof 4+k sits on local edge k, and local edge k connects
local vertices LOCAL_EDGES[k].

Correctness checks, both O(nnz) so they scale to a genuinely large mesh
(unlike building an ndofs_global x ndofs_global dense reference matrix):
  1. Exact match against an independent, basix-quadrature-based
     reference (mirrors test_emit.py's _reference_stiffness, kept
     duplicated rather than imported so demo/ stays free of a test/-
     directory dependency -- see demo/README.md), run on the smallest
     possible mesh (n=1, one cube, 6 cells) where a dense comparison is
     still cheap. This is the first real test that dof-sharing across
     cells (this cube's cells already share several interior
     faces/edges) combines correctly.
  2. The "patch test", on the FULL mesh at whatever -n was given: a
     stiffness matrix applied to the constant function is zero (since
     grad(1) = 0), and the nodal Lagrange basis partitions unity, so
     row i's entries must sum to (numerically) zero for every i. This
     independently validates the entire assembly -- dof numbering,
     per-cell geometry, and the kernel's own arithmetic -- and stays
     cheap (one pass over nnz) no matter how large the mesh gets.
  3. Symmetry: A[i, j] == A[j, i] for every stored entry (also O(nnz)).

Why not the real GPU path (generate_csr_entry_gpu_module): its own test
(test_gpu_assembly.py's
test_generate_csr_entry_gpu_module_matches_quadrature_reference_via_execution_engine)
needs a real NVIDIA GPU and a CUDA-enabled MLIR build, and explicitly
skips (rather than failing) without one -- "e.g. on the Mac this was
developed on, which has no CUDA Toolkit at all". That path is also
currently scoped to one cell per gpu.launch_func call (see that
function's docstring) rather than a whole mesh per launch, so even with
CUDA available it would still mean one Python-level call per cell --
multi-cell batching in a single launch is explicitly future work there.
This script's CPU path, by contrast, already needs only the one call
built by generate_csr_assembly_module below; the GPU path's remaining
advantage over it is running the per-cell work in hardware parallel, not
reducing how many times Python has to call into the kernel.

Usage:
    python3 demo/assemble_mesh_gpu.py [degree] [n]

    degree: Lagrange degree, default 2. Only 1 and 2 are supported --
        generate_csr_assembly_module only extracts geometry for affine
        tetrahedra (any degree), but the dof map built here only knows
        how to place vertex dofs (P1) and vertex + edge-midpoint dofs
        (P2); degree 3+ would additionally need face/interior dof
        placement, not implemented in this script.
    n: mesh resolution, default 6 -- n x n x n cubes, 6*n**3 tetrahedra
        (n=6 -> 1296 cells, 2197 dofs, ~55k nonzeros at P2). Assembly
        itself is a SINGLE call into the generated kernel regardless of
        n (see assemble_global_matrix and generate_csr_assembly_module's
        own docstring): runtime is dominated by that kernel's internal
        cell loop actually executing 6*n**3 cell iterations, each doing
        ndofs**2 quadrature-weighted accumulations plus one
        binary-search-and-atomic-add CSR scatter per (i, j) pair, not by
        any Python-level per-call overhead -- there is no longer a
        per-triple Python call to have overhead in the first place. This
        script demonstrates assembly *correctness at real mesh scale*,
        not GPU throughput -- generate_csr_entry_gpu_module (still one
        gpu.launch_func call per cell, see its own docstring) is what
        would move this cell loop onto actual GPU hardware, which this
        script can't exercise here (see above).
"""

from __future__ import annotations

import ctypes
import sys
import time
from pathlib import Path

import basix
import numpy as np
from basix_uflx import element
from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner

from uflx_mlir.gpu_assembly import generate_csr_assembly_module

sys.path.insert(0, str(Path(__file__).parent))
import harness as mlir_harness  # noqa: E402  (see the sys.path.insert above, matches every other demo/*.py script)

CELL = basix.CellType.tetrahedron

# Local edge k (of a tetrahedron's 6 edges) connects local vertices
# LOCAL_EDGES[k]; local dof 4+k sits on that edge for a P2-equispaced
# Lagrange element. Confirmed directly against basix (not guessed):
#   basix.cell.sub_entity_connectivity(basix.CellType.tetrahedron)[1]
#     -> edge 0=(2,3), 1=(1,3), 2=(1,2), 3=(0,3), 4=(0,2), 5=(0,1)
#   basix.create_element(P, tetrahedron, 2, equispaced).entity_dofs
#     -> dim-1 (edge) entities carry local dofs [4],[5],[6],[7],[8],[9]
#        in that same edge order.
LOCAL_EDGES = [(2, 3), (1, 3), (1, 2), (0, 3), (0, 2), (0, 1)]


def _stiffness_form(degree: int):
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    u = TrialFunction(space)
    v = TestFunction(space)
    return inner(grad(u), grad(v)) * dx, e.dim


def _geometry_from_coords(coords: np.ndarray) -> np.ndarray:
    """Packed affine tetrahedral Poisson metric fed to the CSR-entry
    kernel's `geometry` argument -- mirrors
    uflx_mlir.geometry.extract_affine_poisson_geometry / test_emit.py's
    _reference_geometry exactly (upper triangle, row-major, of
    |detJ| * Jinv @ Jinv.T). Kept duplicated here rather than imported,
    same reasoning as _reference_stiffness_cell below."""
    x0, x1, x2, x3 = coords
    jacobian = np.column_stack([x1 - x0, x2 - x0, x3 - x0])
    jacobian_inv = np.linalg.inv(jacobian)
    metric = abs(np.linalg.det(jacobian)) * jacobian_inv @ jacobian_inv.T
    return metric[np.triu_indices(3)]


def _reference_stiffness_cell(coords: np.ndarray, degree: int) -> np.ndarray:
    """Independent basix-quadrature reference local stiffness matrix --
    mirrors test_emit.py's _reference_stiffness exactly. Kept duplicated
    (not imported from test/) so demo/ stays free of a test/-directory
    dependency, matching every other script in this folder."""
    e = basix.create_element(basix.ElementFamily.P, CELL, degree, basix.LagrangeVariant.equispaced)
    qdeg = max(2 * (degree - 1), 1)
    points, weights = basix.make_quadrature(CELL, qdeg)
    points = np.asarray(points, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    tab = np.asarray(e.tabulate(1, points))

    x0, x1, x2, x3 = coords
    jacobian = np.column_stack([x1 - x0, x2 - x0, x3 - x0])
    det_j = np.linalg.det(jacobian)
    jacobian_inv = np.linalg.inv(jacobian)

    ndofs = tab.shape[2]
    local = np.zeros((ndofs, ndofs))
    for q, w in enumerate(weights):
        grads_ref = tab[1:4, q, :, 0].T
        grads_phys = grads_ref @ jacobian_inv
        local += (grads_phys @ grads_phys.T) * abs(det_j) * w
    return local


def _kuhn_triangulate_cube():
    """The 6 tets making up a unit cube, standard Kuhn/Freudenthal
    decomposition: all 6 share the (0,0,0)-(1,1,1) main diagonal, one per
    permutation of the 3 axes. Each tet's own local vertex order (v000,
    ..., v111) becomes local dofs 0..3 in exactly that order once mapped
    to global vertex ids in build_mesh() below -- must stay consistent
    with _geometry_from_coords/_reference_stiffness_cell, which read
    `coords` as x0..x3 in this same order. Orientation (i.e. whether
    (x1-x0, x2-x0, x3-x0) is a right- or left-handed frame) is
    deliberately not controlled for -- see this module's docstring."""
    v000 = (0, 0, 0)
    axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    tets = []
    for perm in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
        corner = list(v000)
        path = [tuple(corner)]
        for axis_idx in perm:
            dx, dy, dz = axes[axis_idx]
            corner = [corner[0] + dx, corner[1] + dy, corner[2] + dz]
            path.append(tuple(corner))
        tets.append(tuple(path))
    return tets


def build_mesh(n: int) -> tuple[np.ndarray, list[tuple[int, int, int, int]]]:
    """A structured mesh of the unit cube: n x n x n sub-cubes, 6 tets each.

    Returns:
        coords: (nverts, 3) float64 vertex coordinates.
        cells: length-6*n**3 list of 4-tuples of global vertex ids.
    """
    npts = n + 1

    def vid(i: int, j: int, k: int) -> int:
        return i * npts * npts + j * npts + k

    coords = np.array(
        [[i / n, j / n, k / n] for i in range(npts) for j in range(npts) for k in range(npts)],
        dtype=np.float64,
    )

    cells = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for tet in _kuhn_triangulate_cube():
                    cells.append(tuple(vid(i + dx, j + dy, k + dz) for dx, dy, dz in tet))
    return coords, cells


def build_dofmap(
    cells: list[tuple[int, int, int, int]], nverts: int, degree: int
) -> tuple[np.ndarray, int, int]:
    """Map each cell's local dofs (basix's own canonical order) to global dof ids.

    P1: local dof i (0..3) is global vertex i of the cell -- no extra dofs.
    P2: local dofs 0..3 as above, plus local dof 4+k on local edge
        LOCAL_EDGES[k] (k=0..5) -- a fresh global id assigned the first
        time each (unordered) global-vertex-pair is seen anywhere in the
        mesh, so cells sharing an edge share its dof id, numbered
        starting right after the vertex dofs.

    Returns:
        (cell_dofs, ndofs, ndofs_global): cell_dofs is a flat int32 array
        (cell*ndofs + local -> global dof id); ndofs is local dofs per
        cell; ndofs_global is the total dof count.

    Raises:
        NotImplementedError: for any degree other than 1 or 2.
    """
    if degree == 1:
        cell_dofs = np.array(cells, dtype=np.int32).reshape(-1)
        return cell_dofs, 4, nverts

    if degree != 2:
        raise NotImplementedError(
            f"build_dofmap only supports degree 1 or 2 (got {degree}) -- "
            "degree 3+ needs face/interior dof placement this mesh generator "
            "does not implement."
        )

    edge_id: dict[tuple[int, int], int] = {}
    cell_dofs = np.empty((len(cells), 10), dtype=np.int32)
    for c, verts in enumerate(cells):
        cell_dofs[c, 0:4] = verts
        for local_k, (a, b) in enumerate(LOCAL_EDGES):
            key = (verts[a], verts[b]) if verts[a] < verts[b] else (verts[b], verts[a])
            gid = edge_id.get(key)
            if gid is None:
                gid = nverts + len(edge_id)
                edge_id[key] = gid
            cell_dofs[c, 4 + local_k] = gid
    ndofs_global = nverts + len(edge_id)
    return cell_dofs.reshape(-1), 10, ndofs_global


def build_csr_pattern(
    cell_dofs: np.ndarray, ndofs: int, ncells: int, ndofs_global: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The global sparsity pattern: row r's columns are the union, over
    every cell touching global dof r, of that cell's other local dofs
    (including r itself, for the diagonal).

    Columns are sorted ascending within each row -- REQUIRED, not just
    tidy: both generate_csr_entry_module's and generate_csr_assembly_module's
    kernels do their own binary search over each row's stored columns (see
    gpu_assembly.py's _binary_search_and_scatter, shared by both), which
    only works on sorted input.

    Returns:
        (avals, acols, arowptr): avals is zero-initialized float64 (one
        entry per stored (row, col) pair), acols/arowptr are int32.
    """
    rows: list[set[int]] = [set() for _ in range(ndofs_global)]
    cell_dofs2d = cell_dofs.reshape(ncells, ndofs)
    for verts in cell_dofs2d:
        verts_list = [int(v) for v in verts]
        for a in verts_list:
            rows[a].update(verts_list)

    arowptr = np.zeros(ndofs_global + 1, dtype=np.int32)
    acols_list: list[int] = []
    for r, cols in enumerate(rows):
        sorted_cols = sorted(cols)
        acols_list.extend(sorted_cols)
        arowptr[r + 1] = arowptr[r] + len(sorted_cols)
    acols = np.array(acols_list, dtype=np.int32)
    avals = np.zeros(len(acols_list), dtype=np.float64)
    return avals, acols, arowptr


def assemble_global_matrix(
    coords: np.ndarray,
    cells: list[tuple[int, int, int, int]],
    cell_dofs: np.ndarray,
    ndofs: int,
    ndofs_global: int,
    degree: int,
    kernel_name: str = "tabulate_tensor_csr_assemble",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build generate_csr_assembly_module's single whole-mesh kernel for
    `degree` and call it exactly ONCE to assemble the full global CSR
    stiffness matrix for the given mesh.

    Unlike an earlier version of this function (one Python call per
    (cell, tx, ty) triple via generate_csr_entry_module), the kernel built
    here takes the whole mesh's packed geometry and flat dof map as plain
    arguments and loops over every cell itself -- see
    generate_csr_assembly_module's own docstring for how its internal
    cell loop, per-cell geometry refresh, and zero-init-then-accumulate
    buffer reuse work. There is no host-side loop over cells or dof pairs
    left in this function at all.

    Returns:
        (avals, acols, arowptr): the assembled CSR matrix.
    """
    form, ndofs_check = _stiffness_form(degree)
    if ndofs_check != ndofs:
        raise AssertionError(
            f"local dof count mismatch: build_dofmap says {ndofs}, "
            f"the P{degree} element says {ndofs_check}"
        )
    module, layout = generate_csr_assembly_module(form, degree, kernel_name, CELL)
    if layout.ndofs != ndofs:
        raise AssertionError(f"layout.ndofs={layout.ndofs} != ndofs={ndofs}")

    ncells = len(cells)
    avals, acols, arowptr = build_csr_pattern(cell_dofs, ndofs, ncells, ndofs_global)

    # Precompute every cell's packed geometry up front, flattened
    # row-major over cells: cell c's own layout.geometry_size-element
    # block sits at geometries[c*layout.geometry_size:(c+1)*layout.geometry_size]
    # -- exactly the layout generate_csr_assembly_module's `geometries`
    # argument expects (see that function's docstring). cell_dofs is
    # already flat in that same per-cell layout (see build_dofmap).
    from mlir.runtime import get_ranked_memref_descriptor

    geometries = np.empty((ncells, layout.geometry_size), dtype=np.float64)
    for c, verts in enumerate(cells):
        geometries[c] = _geometry_from_coords(coords[list(verts)])
    geometries = geometries.reshape(-1)

    engine = mlir_harness.build_engine_from_module(module, mlir_harness.UFLX_PIPELINE)
    raw_fn = engine.lookup(kernel_name)

    avals_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(avals)))
    acols_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(acols)))
    arowptr_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(arowptr)))
    geometries_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(geometries)))
    cell_dofs_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(cell_dofs)))
    ncells_p = ctypes.pointer(ctypes.c_longlong(ncells))

    packed = (ctypes.c_void_p * 6)(
        ctypes.cast(avals_pp, ctypes.c_void_p).value,
        ctypes.cast(acols_pp, ctypes.c_void_p).value,
        ctypes.cast(arowptr_pp, ctypes.c_void_p).value,
        ctypes.cast(geometries_pp, ctypes.c_void_p).value,
        ctypes.cast(cell_dofs_pp, ctypes.c_void_p).value,
        ctypes.cast(ncells_p, ctypes.c_void_p).value,
    )

    t0 = time.perf_counter()
    raw_fn(packed)  # ONE call assembles the entire mesh -- no Python-level
    # loop over cells or dof pairs; see generate_csr_assembly_module's
    # docstring for what happens inside this single call.
    t1 = time.perf_counter()

    print(
        f"  assembled: {ncells} cells, {ndofs_global} dofs, {len(acols)} nonzeros, "
        f"1 kernel call in {t1 - t0:.3f}s"
    )
    return avals, acols, arowptr


def check_small_mesh_against_reference(degree: int) -> None:
    """Exact check on the smallest possible mesh (1 cube, 6 cells) against
    an independent basix-quadrature reference -- see module docstring,
    check (1)."""
    coords, cells = build_mesh(1)
    cell_dofs, ndofs, ndofs_global = build_dofmap(cells, len(coords), degree)
    avals, acols, arowptr = assemble_global_matrix(
        coords, cells, cell_dofs, ndofs, ndofs_global, degree
    )

    a = np.zeros((ndofs_global, ndofs_global))
    for r in range(ndofs_global):
        for idx in range(arowptr[r], arowptr[r + 1]):
            a[r, acols[idx]] = avals[idx]

    a_ref = np.zeros((ndofs_global, ndofs_global))
    cell_dofs2d = cell_dofs.reshape(len(cells), ndofs)
    for c, verts in enumerate(cells):
        local = _reference_stiffness_cell(coords[list(verts)], degree)
        gdofs = cell_dofs2d[c]
        a_ref[np.ix_(gdofs, gdofs)] += local

    np.testing.assert_allclose(a, a_ref, rtol=1e-9, atol=1e-8)
    print(f"  P{degree} exact check (1 cube, {len(cells)} cells): MATCH")


def patch_test(avals: np.ndarray, acols: np.ndarray, arowptr: np.ndarray, ndofs_global: int) -> None:
    """Row sums must vanish -- see module docstring, check (2). O(nnz)."""
    row_of_entry = np.repeat(np.arange(ndofs_global), np.diff(arowptr))
    row_sums = np.bincount(row_of_entry, weights=avals, minlength=ndofs_global)
    max_abs = float(np.max(np.abs(row_sums)))
    print(f"  patch test (row sums should be ~0): max |row sum| = {max_abs:.3e}")
    assert max_abs < 1e-6, "patch test failed -- assembly is wrong somewhere"


def check_symmetry(avals: np.ndarray, acols: np.ndarray, arowptr: np.ndarray, ndofs_global: int) -> None:
    """A[i, j] == A[j, i] for every stored entry -- see module docstring,
    check (3). O(nnz) via a dict keyed by (row, col)."""
    entries: dict[tuple[int, int], float] = {}
    for r in range(ndofs_global):
        for idx in range(arowptr[r], arowptr[r + 1]):
            entries[(r, int(acols[idx]))] = float(avals[idx])
    max_asym = 0.0
    for (r, c), v in entries.items():
        max_asym = max(max_asym, abs(v - entries.get((c, r), 0.0)))
    print(f"  symmetry check: max |A[i,j] - A[j,i]| = {max_asym:.3e}")
    assert max_asym < 1e-8, "symmetry check failed -- assembly is wrong somewhere"


def main() -> None:
    degree = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 6

    print(f"--- exact correctness check (P{degree}, 1x1x1 cube) ---")
    check_small_mesh_against_reference(degree)

    ncells = 6 * n**3
    print(f"\n--- P{degree} assembly, {n}x{n}x{n} mesh ({ncells} cells) ---")
    coords, cells = build_mesh(n)
    cell_dofs, ndofs, ndofs_global = build_dofmap(cells, len(coords), degree)
    avals, acols, arowptr = assemble_global_matrix(
        coords, cells, cell_dofs, ndofs, ndofs_global, degree
    )
    patch_test(avals, acols, arowptr, ndofs_global)
    check_symmetry(avals, acols, arowptr, ndofs_global)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
