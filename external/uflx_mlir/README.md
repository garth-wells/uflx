# uflx-mlir

An MLIR code-generation backend for UFLx forms: an alternative to
`uflx_codegeneration`'s C backend that builds MLIR IR directly (via the
Python op-builder API) and JITs it in-process with
`mlir.execution_engine.ExecutionEngine`, instead of generating C and
shelling out to a system compiler.

`generate_mlir_module(form, degree, kernel_name, cell)` reuses
`uflx_codegeneration`'s own lowering pipeline (quadrature
selection/tabulation, geometry expansion, pull-back/push-forward,
inner-product expansion), adds the experimental extraction described below,
and directly expands any remaining affine tetrahedral geometry cellwise. It
then replaces the final "walk the lowered graph and emit code" step and
reorders and hoists the generated loop
nest (`uflx_mlir/hoist.py`) -- `uflx_codegeneration`'s pipeline nests the
per-dof loops OUTSIDE the quadrature loop, which is the worst order for
this kind of assembly, since it forces recomputing quadrature-point-only
terms (most importantly the cell's Jacobian/detJ/cofactors) once per
(dof, dof, quadrature-point) combination instead of once per quadrature
point. See `hoist.py`'s module docstring for the numbers: on a P3
tetrahedron stiffness kernel, this took the generated kernel from ~27x
slower than FFCx's own compiled kernel to ~1.8x faster, per call
(measured on one Apple Silicon laptop -- not a general performance claim
across hardware, but a real, JIT'd, end-to-end measurement, not an
estimate).

For affine tetrahedral Poisson forms, the backend also extracts the
cellwise geometry tensor

```text
G = abs(det(J)) * inv(J) * inv(J).T
```

into a separate exported function named `<kernel_name>_geometry`. The
function accepts a `memref<6xf64>` output followed by the usual
`memref<4x3xf64>` coordinate dofs and stores the symmetric tensor in the
order `(G00, G01, G02, G11, G12, G22)`. The tabulation function keeps its
existing public signature, calls this geometry function once per cell,
and consumes the six resulting values in its reference-gradient
contraction. `geometry_kernel_name(kernel_name)` returns the exported
symbol name for clients that want to invoke the geometry kernel directly.
Passing `inline_geometry=True` to `generate_mlir_module` instead computes
the six components directly in the assembly function, avoiding the call and
temporary geometry buffer while retaining the separately callable exported
geometry function in the module.

Fission scratch buffers are stack-allocated once in the function entry
block. When equal test and trial spaces produce expressions that differ
only by their dof-loop variable, the generated kernel reuses the same
scratch values for both sides. This avoids both hot-loop heap allocation
and duplicate basis transformations without assuming that all forms are
Galerkin forms.

## Status / restrictions

Same restrictions as `uflx_codegeneration`'s own pipeline currently has: a
single cell (codim-0) integral, no coefficients or constants, one scalar
element per function space (see `integrals_to_quadrature`'s own asserts).
Geometry-kernel extraction is currently deliberately narrower: it covers
only scalar Poisson stiffness forms on affine tetrahedra. Other forms use
the existing inline geometry lowering unchanged.

## Installing

```bash
pip install -e external/uflx_mlir
```

This installs `uflx-mlir`'s own dependencies (`uflx`, `uflx-codegeneration`,
`numpy`, `fenics-basix`), same as the sibling `external/` packages. It does
**not** install the `mlir` Python package itself: no pip wheel reliably
provides Python bindings that are ABI-compatible with a specific LLVM/MLIR
checkout, so it has to be built locally -- see below.

## Building LLVM/MLIR with Python bindings

Common prerequisites (any platform):

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install numpy pybind11 nanobind PyYAML ninja cmake
```

`pip install pybind11` alone can grab a version newer than the LLVM
checkout's pin and break the bindings build -- if the build below fails
inside `mlir/python`, check `llvm-project/mlir/python/requirements.txt` in
your checkout for the exact pinned version and install that instead.

Clone a released branch rather than `main`, for reproducibility -- this
backend was developed and validated against `release/18.x` specifically:

```bash
git clone --depth 1 --branch release/18.x https://github.com/llvm/llvm-project.git
cd llvm-project && mkdir build && cd build
```

### macOS (Apple Silicon / arm64)

Prerequisites: Xcode command line tools (`xcode-select --install`),
`brew install cmake ninja`.

```bash
cmake -G Ninja ../llvm \
  -DLLVM_ENABLE_PROJECTS="mlir" \
  -DLLVM_TARGETS_TO_BUILD="AArch64" \
  -DCMAKE_BUILD_TYPE=Release \
  -DLLVM_ENABLE_ASSERTIONS=ON \
  -DLLVM_BUILD_LLVM_DYLIB=ON \
  -DLLVM_LINK_LLVM_DYLIB=ON \
  -DMLIR_ENABLE_BINDINGS_PYTHON=ON \
  -DPython3_EXECUTABLE="$(which python3)" \
  -DLLVM_INCLUDE_TESTS=OFF \
  -DLLVM_INCLUDE_EXAMPLES=OFF \
  -DLLVM_INCLUDE_BENCHMARKS=OFF \
  -DLLVM_OPTIMIZED_TABLEGEN=ON

ninja -j4   # keep modest on unified-memory Macs -- a full -j$(nproc) can OOM while linking
```

### macOS (Intel / x86_64)

Identical to the above with `-DLLVM_TARGETS_TO_BUILD="X86"`. `-j` can
usually go higher than 4 (check `sysctl -n hw.ncpu`) since Intel Macs
don't share Apple Silicon's unified memory pool as tightly -- still watch
memory usage during the link step.

### Linux (Ubuntu/Debian, x86_64 or aarch64)

Prerequisites:

```bash
sudo apt-get install build-essential cmake ninja-build python3-dev python3-venv
sudo apt-get install clang lld   # optional, faster build/link than the default toolchain
```

```bash
cmake -G Ninja ../llvm \
  -DLLVM_ENABLE_PROJECTS="mlir" \
  -DLLVM_TARGETS_TO_BUILD="host" \
  -DCMAKE_BUILD_TYPE=Release \
  -DLLVM_ENABLE_ASSERTIONS=ON \
  -DLLVM_BUILD_LLVM_DYLIB=ON \
  -DLLVM_LINK_LLVM_DYLIB=ON \
  -DMLIR_ENABLE_BINDINGS_PYTHON=ON \
  -DPython3_EXECUTABLE="$(which python3)" \
  -DLLVM_INCLUDE_TESTS=OFF \
  -DLLVM_INCLUDE_EXAMPLES=OFF \
  -DLLVM_INCLUDE_BENCHMARKS=OFF \
  -DLLVM_OPTIMIZED_TABLEGEN=ON \
  -DLLVM_USE_LINKER=lld

ninja -j$(nproc)
```

`-DLLVM_TARGETS_TO_BUILD="host"` auto-detects the machine's own
architecture (x86_64 or aarch64) rather than hardcoding one -- use an
explicit list like `"X86;AArch64"` instead if the built compiler needs to
target more than just the machine it's built on. Drop
`-DLLVM_USE_LINKER=lld` if `lld` isn't installed; the default linker
works, just slower.

### After building (any platform)

```bash
export PYTHONPATH=$PWD/tools/mlir/python_packages/mlir_core:$PYTHONPATH
python3 -c "import mlir; print(mlir.__file__)"   # should print a path, no error
```

Add the `export PYTHONPATH=...` line to your shell profile (or a `.env`
you source), since it's needed every session -- `pip install
external/uflx_mlir` does not put `mlir` on your path.

### Build time and resources

A Release build of just LLVM+MLIR for one target architecture (this
package doesn't need any of LLVM's other backends) takes roughly 20-40
minutes on a modern laptop; expect longer if targeting multiple
architectures. Linking is the memory-hungry step -- if the build runs out
of memory there, lower ninja's job count (`ninja -j2`) rather than
reaching for a bigger machine first.

## Running the tests

```bash
pip install -e "external/uflx_mlir[ci]"
pytest external/uflx_mlir/test
```

The tests import `mlir` directly and are skipped (not failed) if it isn't
importable -- this is what happens in CI, which does not build LLVM from
source for this experimental backend. Run them locally after completing
the build above to actually exercise the generator end to end.

## Origin

This backend originated as a prototype in a companion repo
(`mlir-kernels`) investigating whether MLIR could replace FFCx's
C-codegen backend for FEM kernel generation. That repo's text-emission
variant of this same generator was never carried over -- only the
op-builder-based generator (`generate_mlir_module`) and its loop-hoisting
analysis (`hoist.py`) were, to keep one supported code path here rather
than two behaviourally-equivalent generators. The rest of `mlir-kernels`
(the hand-written/basix-generated quadrature-loop kernels it started
from, and the FFCx comparison harness built around them) has since been
brought over too, under `demo/` -- see that folder's own README. Two of
those scripts (`demo/ffcx_compare.py`, `demo/ffcx_compare_uflx.py`) need
the optional `fenics-ffcx` dependency, which is why they live in `demo/`
rather than as part of this package's own installable surface.


## GPU linear-form vectors

`generate_linear_assembly_gpu_module` emits batched vector assembly for real
linear forms, including `inner(grad(w), grad(v)) * dx` and `inner(w, v) * dx`
with `w` a runtime `Coefficient`. It reuses the CPU coefficient lowering and
supports multiple coefficients, including different coefficient/test spaces.
The existing CSR matrix generators retain their bilinear-form interface.

```python
from uflx_mlir import generate_linear_assembly_gpu_module
from uflx_mlir.gpu_runtime import assemble_linear_gpu

module, layout = generate_linear_assembly_gpu_module(
    form, degree, "assemble_vector", basix.CellType.tetrahedron
)
seconds = assemble_linear_gpu(
    module,
    layout,
    "assemble_vector",
    coordinates,
    coefficients,
    cell_dofs,
    output,
    backend="cuda",
    chip="sm_89",
)
```

Inputs are C-contiguous NumPy arrays: `coordinates` is float64 with shape
`(ncells, *layout.coordinate_shape)`, `coefficients` is float64 with shape
`(ncells, layout.coefficient_size)`, and `cell_dofs` is int32 with shape
`(ncells, layout.ndofs)`. Within each cell, coefficient DOF blocks are concatenated
in increasing `Coefficient.count` order, exactly as in the CPU emitter. Pass
`(ncells, 0)` coefficients for a coefficient-free form. The writable float64
`output` vector is **accumulated into**, so initialize it to zero for fresh
assembly. The execution helper allocates device memory, transfers inputs,
launches, synchronizes and copies output back. It consumes the MLIR module;
generate a new module before another invocation. The returned time includes only
one action launch and synchronization, with no explicit warm-up. Geometry
preparation, when needed, completes before this timer starts.

Affine tetrahedral stiffness actions use a stored geometry metric by default.
The companion GPU kernel `geometry_kernel_name(name)` computes
`G = abs(det(J)) * inv(J) * inv(J).T` once per cell. Only six float64 values are
stored per cell, in order `(G00, G01, G02, G11, G12, G22)`. Quadrature weights
remain in the action kernel, so the metric is independent of polynomial degree
and quadrature rule. Recompute it when the mesh coordinates change.

For these kernels `layout.geometry_size == 6`, and the action's second device
buffer is the flattened metric instead of coordinates. The companion geometry
kernel takes two flattened float64 memrefs `(metric, coordinates)`; launch it
with one x thread per cell, block `(128, 1, 1)` and grid
`((ncells + 127) // 128, 1, 1)`. Applications managing persistent device buffers
can run geometry setup once and then reuse its output for repeated actions.
`assemble_linear_gpu` handles this setup automatically, or accepts a precomputed
host metric array through `geometry=metric` to skip the setup kernel.

Set `precompute_geometry=False` in the generator (or `--inline-geometry` in the
demo) to retain the original path for comparison. Forms whose geometry cannot
be fully represented by the extracted metric retain their coordinate buffer
and have `layout.geometry_size == 0`.

The default stiffness action now follows a cooperative quadrature schedule.
Each block stages cell coefficients and basis tables in aligned shared memory.
One thread per quadrature point evaluates the coefficient gradient and applies
the metric. After a block-wide barrier, one thread per test DOF contracts
those shared fluxes against the test gradient and makes one global atomic add.
The fixed-size coefficient and quadrature contractions are unrolled. This avoids
recomputing the coefficient gradient for every test DOF.

For this path, block x selects cells and block y selects quadrature points or
DOFs. Threads per cell are `max(ndofs, nquadrature)`. The default P1/P2/P3/P4
blocks are `(64, 4, 1)`, `(25, 10, 1)`, `(12, 20, 1)`, and `(7, 35, 1)`.
All threads participate in the barriers, including padding cells. Override
`cells_per_block` or `target_block_size` to tune the grouping, within the
256-thread limit. Duplicate basis tables share one buffer. The generator falls
back to serial quadrature if its shared-memory estimate exceeds 48 KiB.

Set `cooperative=False` in the generator, or `--serial-quadrature` in the demo,
to compare the previous stored-metric kernel. That path keeps one thread per
test DOF, a serial quadrature loop, and cell grouping in z. When unspecified,
`target_block_size` is 256 for cooperative quadrature and 128 for the serial
path. Other linear forms retain their existing fallback behavior.

The generator currently requires one integral, one test-DOF axis and at most
256 local test DOFs, with the CPU lowering's existing geometry/scalar limitations.
Gradient-action tests cover affine tetrahedra P1 through P4 on CUDA and HIP,
including multiple cells, distinct cell geometries and coefficients, and shared
DOF scattering. The demo mesh builder supports P1 through P4, including orientation-correct
sharing of edge and face nodes and P4 cell-interior nodes:

```bash
python demo/assemble_linear_gpu.py --degree 1 --n 20 --backend cuda --chip sm_89
python demo/assemble_linear_gpu.py --degree 2 --n 20 --backend amd --chip gfx1100
```

Run hardware tests with `UFLX_GPU_BACKEND=cuda` or `UFLX_GPU_BACKEND=amd`, and
optionally `UFLX_GPU_CHIP`, using `pytest test/test_gpu_linear.py`. Hardware
execution tests skip unless a backend is explicitly selected; generation and
input-validation tests run wherever MLIR Python bindings are available.

### HIP runtime isolation

Both `demo/assemble_mesh_gpu.py` and `assemble_linear_gpu` execute HIP kernels
in a fresh Python worker. MLIR compilation remains in the parent; only HSACO
bytes, launch metadata and NumPy arrays cross the process boundary. The worker
imports neither MLIR nor UFLx. This avoids duplicate LLVM option registration
when the ROCm runtime/COMGR and the MLIR bindings use conflicting LLVM libraries.
No `RTLD_DEEPBIND` or library preloading is needed.

Arrays are exchanged through temporary `.npy` files; nonzero initial output is
preserved and updated output is copied back only after successful completion.
Worker startup and file exchange increase total demo wall time but are excluded
from reported kernel timing, as are device allocation, transfers and geometry
preparation. Worker failures are returned as Python exceptions with stderr.
The worker uses the current Python interpreter and requires NumPy installed in
that environment; it ignores `PYTHONPATH` via Python's isolated mode.
