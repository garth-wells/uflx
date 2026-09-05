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
kernels/                 # generated *.mlir files land here (gitignored)
```

`harness.py`, `generate_kernel.py`, `check_lowering.py`, `run_higher_order.py`,
and `run_uflx_stiffness.py` need only `uflx_mlir`'s own dependencies (basix,
numpy, the MLIR Python bindings) -- see the package README's "Installing"
and "Building LLVM/MLIR with Python bindings" sections. `ffcx_compare.py`
and `ffcx_compare_uflx.py` additionally need:

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
  both, and a programmatic (no shelled-out `mlir-opt` binary) parse+lower
  check.
- Doesn't cover yet: coefficients/constants arguments, facet integrals,
  multiple cells batched in one call, or non-affine (curved) geometry.
- `generate_kernel.py`'s hand-written kernels assume positively-oriented
  tetrahedra (use `detJ` directly rather than `abs(detJ)`) -- fine for a
  hand-built test cell, not for arbitrary DOLFINx meshes, which can have
  either orientation. The UFLx-driven path (`run_uflx_stiffness.py`,
  `ffcx_compare_uflx.py`) does not have this restriction: UFLx's
  `JacobianDeterminant` always takes `abs()` to handle either orientation.
- Degree > 6 equispaced Lagrange nodes get increasingly ill-conditioned
  (Runge-phenomenon-style); `generate_kernel.py` prints a note at that
  point but doesn't stop you.
