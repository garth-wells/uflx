# demo

Worked examples and comparisons for `uflx_mlir`'s GPU/CPU codegen, migrated
in from a companion prototype repo (`mlir-kernels`) once that repo's
generator was superseded by `../uflx_mlir/emit.py`'s own op-builder path
(see the package README's "Origin" section). Kept as a separate `demo/`
folder, rather than merged into the package proper, because two of these
scripts need the optional `fenics-ffcx` dependency, which the installable
`uflx_mlir` package itself does not require.

Cell type is tetrahedron (3D) throughout.

## Layout

```
harness.py              # core JIT harness: parse/lower/JIT/call any kernel, numpy P1 reference
generate_kernel.py       # basix-driven generator for arbitrary-degree (P1, P2, P3, ...) quadrature-loop kernels
check_lowering.py        # pure parse+lower sanity check via the Python bindings, no ExecutionEngine/JIT, no mlir-opt binary needed
run_higher_order.py      # validates a generated P{N} kernel (from generate_kernel.py) against its own reference
run_uflx_stiffness.py    # UFLx form -> ../uflx_mlir/emit.py's op-builder API -> JIT, for any Lagrange degree
ffcx_compare.py          # FFCx-vs-MLIR compile/call-time comparison for a basix-generated kernel (needs fenics-ffcx)
ffcx_compare_uflx.py     # full-pipeline FFCx-vs-UFLx->MLIR comparison (needs fenics-ffcx)
assemble_mesh_gpu.py     # assembles a real global CSR matrix (uflx_mlir.gpu_assembly) on a structured tet mesh, P1/P2
kernels/                 # generated *.mlir files land here (gitignored)
```

`harness.py`, `generate_kernel.py`, `check_lowering.py`, `run_higher_order.py`,
`run_uflx_stiffness.py`, and `assemble_mesh_gpu.py` need only `uflx_mlir`'s
own dependencies (basix, numpy, the MLIR Python bindings) -- see the
package README's "Installing" and "Building LLVM/MLIR with Python
bindings" sections. `ffcx_compare.py` and `ffcx_compare_uflx.py`
additionally need:

```bash
pip install fenics-ffcx
```

## Run

```bash
python3 demo/generate_kernel.py 1        # generate kernels/p1_stiffness.mlir
python3 demo/harness.py                  # P1 kernel: JIT, validate, time
python3 demo/generate_kernel.py 3        # generate kernels/p3_stiffness.mlir
python3 demo/check_lowering.py 3         # sanity-check it parses/lowers (no JIT)
python3 demo/run_higher_order.py 3       # JIT it, validate against basix-derived reference

python3 demo/run_uflx_stiffness.py 3     # UFLx form -> MLIR (op-builder API) -> JIT, degree 3

python3 demo/assemble_mesh_gpu.py 2 6     # assemble a real P2 global CSR matrix, 6x6x6 mesh (1296 cells)

# needs fenics-ffcx (see above):
python3 demo/ffcx_compare.py 3           # compile-time + call-time vs real FFCx codegen
python3 demo/ffcx_compare_uflx.py 3      # same, but for the full UFLx->MLIR pipeline
```

`ffcx_compare.py`'s default (no arg) is degree 1. For degree > 1, run
`generate_kernel.py <degree>` first so the kernel file exists.
`run_uflx_stiffness.py` and `ffcx_compare_uflx.py` build their kernel
straight from a UFLx form via `../uflx_mlir/emit.py`, so they need no
separate generation step.

## Results so far

- P1 (generated one-point quadrature loop): MLIR's direct-ctypes calling
  convention beats FFCx per-call, and MLIR's JIT compile is much faster
  than FFCx's C compile -- both compile time and cold-start total favor
  MLIR by roughly an order of magnitude at this degree.
- P2/P3 (quadrature-loop kernels): still favor MLIR, once the FFCx-side
  Lagrange variant is matched (`equispaced`, matching `LAGRANGE_VARIANT`
  in `generate_kernel.py` -- degree <=2 has <=1 point per edge so any
  variant coincides, but degree >=3 doesn't, and a variant mismatch shows
  up as a silent wrong-answer, not a crash).
  - The naive generated kernel recomputed per-quadrature-point gradient
    terms redundantly inside the dof-pair double loop; precomputing them
    once per quadrature point into scratch buffers (see
    `generate_kernel.py`'s `numx_scratch`/`numy_scratch`/`numz_scratch` and
    the precompute `scf.for` loop) removed that redundant work.
- P4/P5: before the precompute-scratch fix, FFCx was noticeably faster
  (roughly 2-3x) per call at these higher degrees, where the naive
  redundant-recompute cost scales up fastest.
- UFLx form -> MLIR (`ffcx_compare_uflx.py`): the naive UFLx-driven P3
  stiffness kernel initially ran ~27x slower per call than FFCx, because
  the whole per-entry expression -- including the Jacobian/detJ/cofactor
  geometry terms, invariant across every (i, j, quadrature-point)
  combination -- was flattened into the innermost loop body with no
  per-node reuse across loop levels. `../uflx_mlir/emit.py`'s
  `generate_mlir_module` now does the same loop-hoisting/fission
  restructuring `generate_kernel.py`'s hand-written kernels always did
  (see `../uflx_mlir/hoist.py`), closing most of that gap; the `[licm]`
  pipeline variant in `ffcx_compare_uflx.py`'s output is kept as a sanity
  check that MLIR's own loop-invariant-code-motion doesn't already do
  `hoist.py`'s job for free, and the `[baseline]` row is kept for
  historical comparison.

## Calling convention notes

`ExecutionEngine.invoke()` works but has real per-call overhead (symbol
lookup + argument marshaling each time). `ExecutionEngine.lookup(name)`
returns a bound ctypes callable using MLIR's generic "packed" convention:
one argument, a pointer to an array of `void*`, where slot *i* holds the
**address of** the i-th real argument -- not the argument's value, and
not one ctypes arg per kernel argument. Getting this wrong is a
segfault, not a clean error; `harness.build_caller_for()` is the
confirmed-correct, reusable implementation. Prefer it -- `run_kernel()`
and `ffcx_compare.py`'s `direct_ctypes` variant both use it by default.

## What this does and doesn't cover yet

- Covers: generate -> lower -> JIT -> call for generalized basix-driven
  quadrature-loop kernels (P1 and up, arbitrary degree), a UFLx-form ->
  MLIR path via `../uflx_mlir/emit.py`, real FFCx-codegen comparisons for
  both, a programmatic (no shelled-out `mlir-opt` binary) parse+lower
  check, and (`assemble_mesh_gpu.py`) assembling a real global CSR matrix
  with a SINGLE call into `../uflx_mlir/gpu_assembly.py`'s whole-mesh
  assembly kernel (`generate_csr_assembly_module`, which loops over every
  cell itself) on a genuine (if modest-sized) P1/P2 tetrahedral mesh with
  real shared dofs -- checked against an independent quadrature reference
  on the smallest mesh, and via a patch-test (row sums vanish) plus a
  symmetry check at whatever mesh size is asked for.
- Doesn't cover yet: coefficients/constants arguments, facet integrals,
  multiple cells batched in one *GPU-launch* call (the CPU-callable
  `generate_csr_assembly_module` path already batches every cell into one
  call -- see above), non-affine (curved) geometry, or
  (`assemble_mesh_gpu.py` specifically) degree > 2 dof placement (no
  face/interior dofs) or actually launching the GPU-wrapped kernel
  (`generate_csr_entry_gpu_module`) on real hardware -- that needs a
  CUDA-capable machine this package wasn't developed against; see
  `assemble_mesh_gpu.py`'s own docstring.
- `generate_kernel.py`'s hand-written kernels assume positively-oriented
  tetrahedra (use `detJ` directly rather than `abs(detJ)`) -- fine for a
  hand-built test cell, not for arbitrary DOLFINx meshes, which can have
  either orientation. The UFLx-driven path (`run_uflx_stiffness.py`,
  `ffcx_compare_uflx.py`) does not have this restriction: UFLx's
  `JacobianDeterminant` always takes `abs()` to handle either orientation.
- Degree > 6 equispaced Lagrange nodes get increasingly ill-conditioned
  (Runge-phenomenon-style); `generate_kernel.py` prints a note at that
  point but doesn't stop you.
