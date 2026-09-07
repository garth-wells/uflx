"""Assemble a global CSR stiffness matrix using CPU, CUDA, or AMD.

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

GPU backends: this script can also assemble via the actual
GPU-launched kernel (generate_csr_assembly_gpu_module) instead of
generate_csr_assembly_module's CPU path -- pass "cuda" (or the legacy
alias "gpu") or "amd" as this script's third argument (see Usage below).
The CUDA path needs a real NVIDIA GPU and an MLIR build configured with
-DMLIR_ENABLE_CUDA_RUNNER=ON (see
test_generate_csr_assembly_gpu_module_matches_quadrature_reference_two_cells's
own docstring in test/test_gpu_assembly.py for exactly what that means
and how libmlir_cuda_runtime.so is located) -- not available on the Mac
this script was originally developed on, which has no CUDA Toolkit at
all. That kernel itself is confirmed correct on a real CUDA machine
(eng-nvidia) at genuine two-cell scale (see test_gpu_assembly.py's own
test_generate_csr_assembly_gpu_module_matches_quadrature_reference_two_cells);
this script's own full-mesh GPU run is what actually exercises it at
real mesh scale. Unlike generate_csr_entry_gpu_module
(the older, one-cell-per-gpu.launch_func-call kernel this module also
still provides -- see its own docstring; deliberately NOT wired into
this script, since a whole mesh's worth of launches that way would be
far slower than either backend this script does use), the batched kernel
used here needs only ONE gpu.launch_func call regardless of mesh size:
gridDim.x = ncells (one block per cell), blockDim = (ndofs, ndofs, 1)
(one thread per local (i, j) entry, same as the older kernel) -- no
Python-level loop over cells on the GPU backend either, matching the CPU
backend's own single-call design.

The AMD path lowers the same gpu.module through ROCDL, asks ROCm's own
clang and ld.lld to build an HSACO code object, and launches it through
the HIP module API. It auto-detects the GPU architecture with
``/opt/rocm/bin/offload-arch`` unless a target chip is supplied.

Usage:
    python3 demo/assemble_mesh_gpu.py [degree] [n] [backend] [target_chip]

    degree: Lagrange degree, default 2. Only 1 and 2 are supported --
        generate_csr_assembly_module/generate_csr_assembly_gpu_module only
        extract geometry for affine tetrahedra (any degree), but the dof
        map built here only knows how to place vertex dofs (P1) and
        vertex + edge-midpoint dofs (P2); degree 3+ would additionally
        need face/interior dof placement, not implemented in this
        script.
    n: mesh resolution, default 6 -- n x n x n cubes, 6*n**3 tetrahedra
        (n=6 -> 1296 cells, 2197 dofs, ~55k nonzeros at P2). Assembly
        itself is a SINGLE kernel call/launch regardless of n or
        backend (see assemble_global_matrix/assemble_global_matrix_gpu
        and generate_csr_assembly_module's/generate_csr_assembly_gpu_module's
        own docstrings): on the CPU backend, runtime is dominated by
        that kernel's internal cell loop actually executing 6*n**3 cell
        iterations; on the GPU backend, those same 6*n**3 cells instead
        run as gridDim.x blocks in real hardware parallel. Either way
        there is no Python-level per-cell or per-triple call overhead --
        this script demonstrates assembly *correctness at real mesh
        scale* on whichever backend is asked for, not a throughput
        comparison between them.
    backend: "cpu" (default), "cuda"/"gpu", or "amd" -- see "GPU
        backends" above. A requested accelerator backend raises a clear
        error rather than silently falling back to CPU.
    target_chip: optional target architecture. CUDA defaults to "sm_80";
        AMD auto-detects it with offload-arch. Examples are "sm_89" for
        eng-nvidia and "gfx1100" for eng-amd. Ignored by the CPU backend.
"""

from __future__ import annotations

import ctypes
import glob
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import basix
import numpy as np
from basix_uflx import element
from uflx import (
    TestFunction,
    TrialFunction,
    coordinate_element,
    dx,
    function_space,
    grad,
    inner,
)

from uflx_mlir.gpu_assembly import (
    assemble_amdgcn_to_hsaco,
    extract_amdgcn_text,
    generate_csr_assembly_gpu_module,
    generate_csr_assembly_module,
    gpu_launch_name,
    lower_module_to_nvvm,
    lower_module_to_rocdl,
)

sys.path.insert(0, str(Path(__file__).parent))
import harness as mlir_harness

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


def _dofs_per_sec(ndofs_global: int, elapsed_seconds: float) -> float:
    """Compute assembly throughput in degrees of freedom per second.

    This is the throughput this module's assembly functions report
    alongside their own wall-clock timing.
    Floors elapsed_seconds at a tiny epsilon rather than risking
    ZeroDivisionError on an implausibly-fast (sub-nanosecond) timer
    reading, which a tiny mesh's single call could in principle hit.
    """
    return ndofs_global / max(elapsed_seconds, 1e-9)


def _stiffness_form(degree: int):
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    u = TrialFunction(space)
    v = TestFunction(space)
    return inner(grad(u), grad(v)) * dx, e.dim


def _geometry_from_coords(coords: np.ndarray) -> np.ndarray:
    """Pack the affine tetrahedral Poisson metric for the CSR kernel.

    This mirrors
    uflx_mlir.geometry.extract_affine_poisson_geometry / test_emit.py's
    _reference_geometry exactly (upper triangle, row-major, of
    |detJ| * Jinv @ Jinv.T). Kept duplicated here rather than imported,
    same reasoning as _reference_stiffness_cell below.
    """
    x0, x1, x2, x3 = coords
    jacobian = np.column_stack([x1 - x0, x2 - x0, x3 - x0])
    jacobian_inv = np.linalg.inv(jacobian)
    metric = abs(np.linalg.det(jacobian)) * jacobian_inv @ jacobian_inv.T
    return metric[np.triu_indices(3)]


def _reference_stiffness_cell(coords: np.ndarray, degree: int) -> np.ndarray:
    """Compute an independent Basix-quadrature local stiffness matrix.

    This mirrors test_emit.py's _reference_stiffness exactly. Kept duplicated
    (not imported from test/) so demo/ stays free of a test/-directory
    dependency, matching every other script in this folder.
    """
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
    """Construct the six-tetrahedron Kuhn triangulation of a unit cube.

    All six tetrahedra share the (0,0,0)-(1,1,1) main diagonal, one per
    permutation of the 3 axes. Each tet's own local vertex order (v000,
    ..., v111) becomes local dofs 0..3 in exactly that order once mapped
    to global vertex ids in build_mesh() below -- must stay consistent
    with _geometry_from_coords/_reference_stiffness_cell, which read
    `coords` as x0..x3 in this same order. Orientation (i.e. whether
    (x1-x0, x2-x0, x3-x0) is a right- or left-handed frame) is
    deliberately not controlled for -- see this module's docstring.
    """
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

    # Vectorized edge numbering: for every (cell, local edge) pair, the
    # canonicalized (lo, hi) global-vertex-pair packs into one int64 key
    # (lo*nverts + hi -- unique since hi < nverts always), so np.unique's
    # own return_inverse gives each distinct edge exactly one id and maps
    # every occurrence of it (however many cells share it) back to that
    # same id, in one vectorized pass -- no Python-level dict, unlike the
    # per-(cell, edge) dict lookup this replaces. The actual id values
    # differ from that dict's first-seen-order numbering (np.unique
    # returns them in sorted-key order instead), but nothing downstream
    # depends on which specific ids a shared edge gets, only that it gets
    # the SAME one everywhere it's touched -- see this module's own
    # correctness checks, all permutation-invariant in the dof numbering.
    cells_arr = np.asarray(cells, dtype=np.int64)
    ncells = cells_arr.shape[0]
    cell_dofs = np.empty((ncells, 10), dtype=np.int32)
    cell_dofs[:, 0:4] = cells_arr

    local_edges = np.asarray(LOCAL_EDGES, dtype=np.int64)  # shape (6, 2)
    edge_v0 = cells_arr[:, local_edges[:, 0]]  # shape (ncells, 6)
    edge_v1 = cells_arr[:, local_edges[:, 1]]  # shape (ncells, 6)
    lo = np.minimum(edge_v0, edge_v1)
    hi = np.maximum(edge_v0, edge_v1)
    edge_key = (lo * nverts + hi).reshape(-1)

    _, inverse = np.unique(edge_key, return_inverse=True)
    nedges = int(inverse.max()) + 1 if inverse.size else 0
    cell_dofs[:, 4:10] = (nverts + inverse).astype(np.int32).reshape(ncells, 6)

    ndofs_global = nverts + nedges
    return cell_dofs.reshape(-1), 10, ndofs_global


def build_csr_pattern(
    cell_dofs: np.ndarray, ndofs: int, ncells: int, ndofs_global: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the global CSR sparsity pattern.

    Row r's columns are the union, over every cell touching global dof r,
    of that cell's other local dofs
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
    # Vectorized sparsity-pattern build: every cell contributes the full
    # ndofs*ndofs grid of (row, col) pairs among its own local dofs (that
    # cell's row r's columns are exactly its own dof list, for every r it
    # owns -- same rule the old rows[a].update(verts_list) loop encoded,
    # just built via broadcasting instead of a per-cell nested Python
    # loop). Packing each (row, col) pair into one int64 key
    # (row*ndofs_global + col) turns "dedupe the pairs contributed by
    # multiple cells, then sort each row's columns ascending" into "sort
    # the keys, then drop adjacent duplicates": sorting keys ascending
    # sorts by row first (since row dominates the key -- col is always
    # < ndofs_global) and column second within a row, which is exactly
    # the CSR contract _binary_search_and_scatter needs (see this
    # function's own docstring).
    #
    # Deliberately NOT np.unique(keys) here: benchmarked at ~38M keys
    # (n=55, the >1e6-dof scale this was written to handle), plain
    # np.unique took 25-40s on this array size in this numpy build --
    # confirmed via a standalone repro with equivalent-size random int64
    # data, so it's not something about these particular keys. A manual
    # np.sort + boolean-diff dedup needs no argsort (only the deduped
    # sorted keys are needed here, not an inverse mapping) and does the
    # same job in a small fraction of the time -- roughly 15s at the
    # same scale in the same environment, vs. what would likely be
    # minutes-to-hours for the original pure-Python set-based loop.
    #
    # Memory scales with ncells*ndofs**2 raw (row, col) pairs before
    # dedup (int64 keys: 8 bytes each) -- a few hundred MB at the
    # >1e6-dof scale this was written to handle, but worth knowing if
    # pushing to a substantially larger mesh still.
    cell_dofs2d = cell_dofs.reshape(ncells, ndofs).astype(np.int64)
    keys_grid = cell_dofs2d[:, :, None] * ndofs_global + cell_dofs2d[:, None, :]
    keys = keys_grid.reshape(-1)

    sorted_keys = np.sort(keys)
    keep = np.empty(sorted_keys.shape[0], dtype=bool)
    keep[0] = True
    np.not_equal(sorted_keys[1:], sorted_keys[:-1], out=keep[1:])
    unique_keys = sorted_keys[keep]

    rows_u = (unique_keys // ndofs_global).astype(np.int32)
    acols = (unique_keys % ndofs_global).astype(np.int32)
    avals = np.zeros(unique_keys.shape[0], dtype=np.float64)

    arowptr = np.zeros(ndofs_global + 1, dtype=np.int32)
    counts = np.bincount(rows_u, minlength=ndofs_global)
    np.cumsum(counts, out=arowptr[1:])
    return avals, acols, arowptr


def assemble_global_matrix(
    coords: np.ndarray,
    cells: list[tuple[int, int, int, int]],
    cell_dofs: np.ndarray,
    ndofs: int,
    ndofs_global: int,
    degree: int,
    kernel_name: str = "tabulate_tensor_csr_assemble",
    return_timing: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Assemble the full matrix with one whole-mesh CPU kernel call.

    Unlike an earlier version of this function (one Python call per
    (cell, tx, ty) triple via generate_csr_entry_module), the kernel built
    here takes the whole mesh's packed geometry and flat dof map as plain
    arguments and loops over every cell itself -- see
    generate_csr_assembly_module's own docstring for how its internal
    cell loop, per-cell geometry refresh, and zero-init-then-accumulate
    buffer reuse work. There is no host-side loop over cells or dof pairs
    left in this function at all.

    Args:
        coords: Mesh vertex coordinates.
        cells: Tetrahedral cell-to-vertex connectivity.
        cell_dofs: Flattened cell-to-global-dof map.
        ndofs: Number of local element degrees of freedom.
        ndofs_global: Number of global degrees of freedom.
        degree: Lagrange polynomial degree.
        kernel_name: Generated MLIR function name.
        return_timing: when True, also return the single kernel call's
            own wall-clock elapsed seconds as a 4th tuple element --
            exactly the t1 - t0 this function's own print statement
            already reports, just returned numerically instead of only
            printed, for callers (e.g. demo/benchmark_throughput.py)
            that sweep many mesh sizes and need the number, not stdout
            to parse.

    Returns:
        (avals, acols, arowptr): the assembled CSR matrix, or
        (avals, acols, arowptr, elapsed_seconds) if return_timing.
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
        f"1 kernel call in {t1 - t0:.3f}s "
        f"({_dofs_per_sec(ndofs_global, t1 - t0):.3e} dofs/sec)"
    )
    if return_timing:
        return avals, acols, arowptr, t1 - t0
    return avals, acols, arowptr


def _find_cuda_runtime_lib() -> str | None:
    """Locate ``libmlir_cuda_runtime.so``.

    This matches test_gpu_assembly.py's execution-engine tests' search (duplicated
    here rather than imported -- see this module's own "kept duplicated"
    convention for staying free of a test/-directory dependency):
    $MLIR_CUDA_RUNTIME_LIB if set, else searched upward from the mlir
    Python package's own install directory (its usual place:
    <build>/lib/libmlir_cuda_runtime.so, a few levels above
    <build>/tools/mlir/python_packages/mlir_core/mlir/).
    """
    cuda_runtime_lib = os.environ.get("MLIR_CUDA_RUNTIME_LIB")
    if cuda_runtime_lib:
        return cuda_runtime_lib

    import mlir

    start_dirs = []
    for p in getattr(mlir, "__path__", []) or []:
        start_dirs.append(os.path.abspath(p))
    if getattr(mlir, "__file__", None):
        start_dirs.append(os.path.dirname(os.path.abspath(mlir.__file__)))

    for start in start_dirs:
        here = start
        for _ in range(8):
            found = glob.glob(os.path.join(here, "lib", "libmlir_cuda_runtime.so*"))
            if found:
                return found[0]
            parent = os.path.dirname(here)
            if parent == here:
                break
            here = parent
    return None


def assemble_global_matrix_gpu(
    coords: np.ndarray,
    cells: list[tuple[int, int, int, int]],
    cell_dofs: np.ndarray,
    ndofs: int,
    ndofs_global: int,
    degree: int,
    kernel_name: str = "tabulate_tensor_csr_assembly_gpu",
    cubin_chip: str = "sm_80",
    return_timing: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Assemble the full matrix with one batched CUDA launch.

    This makes one real gpu.launch_func call,
    gridDim.x = ncells (one block per cell), blockDim = (ndofs, ndofs, 1)
    -- to assemble the full global CSR stiffness matrix on an actual GPU.

    See assemble_global_matrix's own docstring for what's shared with the
    CPU path (dof-count check, CSR pattern construction, per-cell
    geometry precompute); the differences here are all about actually
    reaching a GPU: locating libmlir_cuda_runtime.so, compiling the
    kernel's gpu.module down to real NVVM/PTX via lower_module_to_nvvm,
    and looking the launch wrapper up by gpu_launch_name(kernel_name)
    (generate_csr_assembly_gpu_module's returned module contains that
    host-side wrapper, not the raw gpu.func symbol itself -- see that
    function's own docstring).

    Needs a real NVIDIA GPU and an MLIR build configured with
    -DMLIR_ENABLE_CUDA_RUNNER=ON -- raises RuntimeError with a clear
    message if libmlir_cuda_runtime.so can't be found, rather than
    silently falling back to the CPU path.

    Args:
        coords: Mesh vertex coordinates.
        cells: Tetrahedral cell-to-vertex connectivity.
        cell_dofs: Flattened cell-to-global-dof map.
        ndofs: Number of local element degrees of freedom.
        ndofs_global: Number of global degrees of freedom.
        degree: Lagrange polynomial degree.
        kernel_name: Generated GPU kernel name.
        cubin_chip: the target NVPTX chip generate_csr_assembly_gpu_module's
            compiled PTX targets -- e.g. "sm_89" for eng-nvidia's Ada
            Lovelace GPU (see lower_module_to_nvvm's own docstring).
            Defaults to "sm_80" (Ampere), matching lower_module_to_nvvm's
            own default; pass whatever matches the actual GPU this runs
            on.
        return_timing: when True, also return the single gpu.launch_func
            call's own wall-clock elapsed seconds as a 4th tuple element
            -- see assemble_global_matrix's own return_timing doc, same
            reasoning.

    Returns:
        (avals, acols, arowptr): the assembled CSR matrix, or
        (avals, acols, arowptr, elapsed_seconds) if return_timing.

    Raises:
        RuntimeError: if libmlir_cuda_runtime.so can't be found (no
            CUDA-enabled MLIR build available).
    """
    form, ndofs_check = _stiffness_form(degree)
    if ndofs_check != ndofs:
        raise AssertionError(
            f"local dof count mismatch: build_dofmap says {ndofs}, "
            f"the P{degree} element says {ndofs_check}"
        )

    cuda_runtime_lib = _find_cuda_runtime_lib()
    if not cuda_runtime_lib or not os.path.exists(cuda_runtime_lib):
        raise RuntimeError(
            "libmlir_cuda_runtime.so not found -- rebuild MLIR with "
            "-DMLIR_ENABLE_CUDA_RUNNER=ON, or set MLIR_CUDA_RUNTIME_LIB "
            "to its path (see test_gpu_assembly.py's own execution-engine "
            "tests for the same check)."
        )

    module, layout = generate_csr_assembly_gpu_module(form, degree, kernel_name, CELL)
    if layout.ndofs != ndofs:
        raise AssertionError(f"layout.ndofs={layout.ndofs} != ndofs={ndofs}")

    lower_module_to_nvvm(module, cubin_chip=cubin_chip)

    ncells = len(cells)
    avals, acols, arowptr = build_csr_pattern(cell_dofs, ndofs, ncells, ndofs_global)

    # Same flat, row-major-over-cells layout as assemble_global_matrix's
    # own geometry precompute -- see that function's docstring.
    geometries = np.empty((ncells, layout.geometry_size), dtype=np.float64)
    for c, verts in enumerate(cells):
        geometries[c] = _geometry_from_coords(coords[list(verts)])
    geometries = geometries.reshape(-1)

    from mlir.execution_engine import ExecutionEngine
    from mlir.runtime import get_ranked_memref_descriptor

    with module.context:
        engine = ExecutionEngine(module, opt_level=3, shared_libs=[cuda_runtime_lib])

    raw_fn = engine.lookup(gpu_launch_name(kernel_name))
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
    raw_fn(packed)  # ONE gpu.launch_func call, gridDim.x=ncells, assembles the whole mesh.
    t1 = time.perf_counter()

    print(
        f"  assembled (GPU, {cubin_chip}): {ncells} cells, {ndofs_global} dofs, "
        f"{len(acols)} nonzeros, 1 launch (gridDim.x={ncells}, "
        f"blockDim=({ndofs},{ndofs},1)) in {t1 - t0:.3f}s "
        f"({_dofs_per_sec(ndofs_global, t1 - t0):.3e} dofs/sec)"
    )
    if return_timing:
        return avals, acols, arowptr, t1 - t0
    return avals, acols, arowptr


def assemble_global_matrix_amd(
    coords: np.ndarray,
    cells: list[tuple[int, int, int, int]],
    cell_dofs: np.ndarray,
    ndofs: int,
    ndofs_global: int,
    degree: int,
    kernel_name: str = "tabulate_tensor_csr_assembly_amd",
    chip: str | None = None,
    rocm_path: str | None = None,
    return_timing: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Assemble the whole mesh in one launch on an AMD GPU through HIP.

    MLIR emits AMDGCN assembly through ROCDL. ROCm's own matching clang
    and ld.lld then produce an HSACO code object, which is loaded and
    launched directly through the HIP module API. This split avoids the
    LLVM-version mismatch that occurs on eng-amd, where the MLIR bindings
    use LLVM 18 but ROCm 7.2's bitcode and linker use LLVM 22.
    """
    form, ndofs_check = _stiffness_form(degree)
    if ndofs_check != ndofs:
        raise AssertionError(
            f"local dof count mismatch: build_dofmap says {ndofs}, "
            f"the P{degree} element says {ndofs_check}"
        )

    rocm = Path(rocm_path or os.environ.get("ROCM_PATH", "/opt/rocm"))
    hip_library = rocm / "lib/libamdhip64.so"
    offload_arch = rocm / "bin/offload-arch"
    required = [
        hip_library,
        offload_arch,
        rocm / "llvm/bin/clang",
        rocm / "llvm/bin/ld.lld",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(f"ROCm installation is incomplete; missing: {', '.join(missing)}")

    if chip is None:
        architectures = subprocess.run(
            [str(offload_arch)], check=True, capture_output=True, text=True
        ).stdout.splitlines()
        if not architectures:
            raise RuntimeError("offload-arch found no AMD GPU")
        chip = architectures[0].split(":", maxsplit=1)[0]

    module, layout = generate_csr_assembly_gpu_module(form, degree, kernel_name, CELL)
    if layout.ndofs != ndofs:
        raise AssertionError(f"layout.ndofs={layout.ndofs} != ndofs={ndofs}")
    lower_module_to_rocdl(module, chip=chip, link_device_libraries=False)
    hsaco = assemble_amdgcn_to_hsaco(extract_amdgcn_text(module), chip=chip, toolkit_path=str(rocm))

    ncells = len(cells)
    avals, acols, arowptr = build_csr_pattern(cell_dofs, ndofs, ncells, ndofs_global)
    geometries = np.empty((ncells, layout.geometry_size), dtype=np.float64)
    for c, verts in enumerate(cells):
        geometries[c] = _geometry_from_coords(coords[list(verts)])
    geometries = geometries.reshape(-1)
    host_arrays = [avals, acols, arowptr, geometries, cell_dofs]

    hip = ctypes.CDLL(str(hip_library))

    def bind(name, restype, *argtypes):
        function = getattr(hip, name)
        function.restype = restype
        function.argtypes = list(argtypes)
        return function

    hip_init = bind("hipInit", ctypes.c_int, ctypes.c_uint)
    hip_set_device = bind("hipSetDevice", ctypes.c_int, ctypes.c_int)
    hip_module_load = bind(
        "hipModuleLoad", ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p
    )
    hip_module_get_function = bind(
        "hipModuleGetFunction",
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.c_char_p,
    )
    hip_malloc = bind("hipMalloc", ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), ctypes.c_size_t)
    hip_memcpy = bind(
        "hipMemcpy",
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.c_int,
    )
    hip_module_launch_kernel = bind(
        "hipModuleLaunchKernel",
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.c_uint,
        ctypes.c_uint,
        ctypes.c_uint,
        ctypes.c_uint,
        ctypes.c_uint,
        ctypes.c_uint,
        ctypes.c_uint,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.POINTER(ctypes.c_void_p),
    )
    hip_device_synchronize = bind("hipDeviceSynchronize", ctypes.c_int)
    hip_free = bind("hipFree", ctypes.c_int, ctypes.c_void_p)
    hip_module_unload = bind("hipModuleUnload", ctypes.c_int, ctypes.c_void_p)
    hip_get_error_string = bind("hipGetErrorString", ctypes.c_char_p, ctypes.c_int)

    def check(code: int, operation: str) -> None:
        if code:
            message = hip_get_error_string(code)
            detail = message.decode() if message else f"HIP error {code}"
            raise RuntimeError(f"{operation}: {detail}")

    check(hip_init(0), "hipInit")
    check(hip_set_device(0), "hipSetDevice")
    hip_module = ctypes.c_void_p()
    hip_function = ctypes.c_void_p()
    allocations: list[ctypes.c_void_p] = []
    elapsed = 0.0
    with tempfile.TemporaryDirectory(prefix="uflx-hip-demo-") as directory:
        hsaco_path = Path(directory) / "kernel.hsaco"
        hsaco_path.write_bytes(hsaco)
        check(
            hip_module_load(ctypes.byref(hip_module), os.fsencode(hsaco_path)),
            "hipModuleLoad",
        )
        check(
            hip_module_get_function(ctypes.byref(hip_function), hip_module, kernel_name.encode()),
            "hipModuleGetFunction",
        )

        arguments = []
        try:
            for array in host_arrays:
                device_pointer = ctypes.c_void_p()
                check(hip_malloc(ctypes.byref(device_pointer), array.nbytes), "hipMalloc")
                allocations.append(device_pointer)
                check(
                    hip_memcpy(
                        device_pointer,
                        ctypes.c_void_p(array.ctypes.data),
                        array.nbytes,
                        1,
                    ),
                    "hipMemcpy host-to-device",
                )
                # Rank-1 MLIR memref ABI: allocated pointer, aligned pointer,
                # offset, size, and stride.
                arguments.extend(
                    [
                        ctypes.c_void_p(device_pointer.value),
                        ctypes.c_void_p(device_pointer.value),
                        ctypes.c_int64(0),
                        ctypes.c_int64(array.size),
                        ctypes.c_int64(1),
                    ]
                )
            kernel_parameters = (ctypes.c_void_p * len(arguments))(
                *(
                    ctypes.cast(ctypes.byref(argument), ctypes.c_void_p).value
                    for argument in arguments
                )
            )
            t0 = time.perf_counter()
            check(
                hip_module_launch_kernel(
                    hip_function,
                    ncells,
                    1,
                    1,
                    ndofs,
                    ndofs,
                    1,
                    0,
                    None,
                    kernel_parameters,
                    None,
                ),
                "hipModuleLaunchKernel",
            )
            check(hip_device_synchronize(), "hipDeviceSynchronize")
            elapsed = time.perf_counter() - t0
            check(
                hip_memcpy(
                    ctypes.c_void_p(avals.ctypes.data),
                    allocations[0],
                    avals.nbytes,
                    2,
                ),
                "hipMemcpy device-to-host",
            )
        finally:
            for device_pointer in allocations:
                check(hip_free(device_pointer), "hipFree")
            if hip_module.value:
                check(hip_module_unload(hip_module), "hipModuleUnload")

    print(
        f"  assembled (AMD, {chip}): {ncells} cells, {ndofs_global} dofs, "
        f"{len(acols)} nonzeros, 1 launch (gridDim.x={ncells}, "
        f"blockDim=({ndofs},{ndofs},1)) in {elapsed:.3f}s "
        f"({_dofs_per_sec(ndofs_global, elapsed):.3e} dofs/sec)"
    )
    if return_timing:
        return avals, acols, arowptr, elapsed
    return avals, acols, arowptr


def check_small_mesh_against_reference(
    degree: int, assemble_fn=assemble_global_matrix, **assemble_kwargs
) -> None:
    """Check the smallest mesh against an independent quadrature reference.

    assemble_fn: assemble_global_matrix (default, CPU) or
        assemble_global_matrix_gpu -- see main() for how the requested
        backend picks which is passed in.
    """
    coords, cells = build_mesh(1)
    cell_dofs, ndofs, ndofs_global = build_dofmap(cells, len(coords), degree)
    avals, acols, arowptr = assemble_fn(
        coords, cells, cell_dofs, ndofs, ndofs_global, degree, **assemble_kwargs
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


def patch_test(
    avals: np.ndarray, acols: np.ndarray, arowptr: np.ndarray, ndofs_global: int
) -> None:
    """Row sums must vanish -- see module docstring, check (2). O(nnz)."""
    row_of_entry = np.repeat(np.arange(ndofs_global), np.diff(arowptr))
    row_sums = np.bincount(row_of_entry, weights=avals, minlength=ndofs_global)
    max_abs = float(np.max(np.abs(row_sums)))
    print(f"  patch test (row sums should be ~0): max |row sum| = {max_abs:.3e}")
    assert max_abs < 1e-6, "patch test failed -- assembly is wrong somewhere"


def check_symmetry(
    avals: np.ndarray, acols: np.ndarray, arowptr: np.ndarray, ndofs_global: int
) -> None:
    """Check matrix symmetry for every stored entry.

    This is O(nnz) and fully vectorized (no per-entry Python loop): pack
    each stored (row, col) into one int64 key (row*ndofs_global + col),
    identical to build_csr_pattern's own packing -- since acols is sorted
    ascending within each row and rows are laid out in increasing order
    (build_csr_pattern's CSR contract), the resulting `keys` array is
    itself already sorted ascending overall, so np.searchsorted can look
    up every entry's mirror (col, row) key directly with no separate
    argsort needed.

    build_csr_pattern's own rows[a].update(verts_list)-style construction
    (see that function's docstring) makes the pattern symmetric by
    construction, so every mirror key is expected to actually be present;
    this still checks that explicitly (rather than assuming it) so a
    future change that broke that invariant would fail loudly here
    instead of silently reading a wrong, unrelated entry.
    """
    row_of_entry = np.repeat(np.arange(ndofs_global, dtype=np.int64), np.diff(arowptr))
    cols64 = acols.astype(np.int64)
    keys = row_of_entry * ndofs_global + cols64
    mirror_keys = cols64 * ndofs_global + row_of_entry

    mirror_idx = np.searchsorted(keys, mirror_keys)
    in_bounds = mirror_idx < keys.shape[0]
    found = np.zeros(mirror_idx.shape[0], dtype=bool)
    found[in_bounds] = keys[mirror_idx[in_bounds]] == mirror_keys[in_bounds]
    if not np.all(found):
        raise AssertionError(
            "sparsity pattern is not symmetric -- found a (row, col) entry "
            "with no (col, row) counterpart; build_csr_pattern should make "
            "this impossible, so this points at a real bug there"
        )

    max_asym = float(np.max(np.abs(avals - avals[mirror_idx])))
    print(f"  symmetry check: max |A[i,j] - A[j,i]| = {max_asym:.3e}")
    assert max_asym < 1e-8, "symmetry check failed -- assembly is wrong somewhere"


def main() -> None:
    """Run the requested mesh-assembly backend and correctness checks."""
    degree = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    backend = sys.argv[3] if len(sys.argv) > 3 else "cpu"
    target_chip = sys.argv[4] if len(sys.argv) > 4 else None

    if backend not in ("cpu", "gpu", "cuda", "amd"):
        raise SystemExit(f"backend must be 'cpu', 'cuda'/'gpu', or 'amd', got {backend!r}")

    if backend in ("gpu", "cuda"):
        cuda_runtime_lib = _find_cuda_runtime_lib()
        if not cuda_runtime_lib or not os.path.exists(cuda_runtime_lib):
            print(
                "GPU backend requested but libmlir_cuda_runtime.so was not "
                "found -- rebuild MLIR with -DMLIR_ENABLE_CUDA_RUNNER=ON, or "
                "set MLIR_CUDA_RUNTIME_LIB to its path. Nothing assembled."
            )
            return
        assemble_fn = assemble_global_matrix_gpu
        assemble_kwargs = {"cubin_chip": target_chip or "sm_80"}
    elif backend == "amd":
        assemble_fn = assemble_global_matrix_amd
        assemble_kwargs = {"chip": target_chip}
    else:
        assemble_fn = assemble_global_matrix
        assemble_kwargs = {}

    print(f"--- exact correctness check (P{degree}, 1x1x1 cube, backend={backend}) ---")
    check_small_mesh_against_reference(degree, assemble_fn, **assemble_kwargs)

    ncells = 6 * n**3
    print(f"\n--- P{degree} assembly, {n}x{n}x{n} mesh ({ncells} cells), backend={backend} ---")
    coords, cells = build_mesh(n)
    cell_dofs, ndofs, ndofs_global = build_dofmap(cells, len(coords), degree)
    avals, acols, arowptr = assemble_fn(
        coords, cells, cell_dofs, ndofs, ndofs_global, degree, **assemble_kwargs
    )
    patch_test(avals, acols, arowptr, ndofs_global)
    check_symmetry(avals, acols, arowptr, ndofs_global)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
