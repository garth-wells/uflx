# Calling an AOT-compiled MLIR kernel from C++

This demo compiles a tetrahedron stiffness kernel to an ordinary object file
and calls it from C++ with two bare pointers:

```cpp
extern "C" void tabulate_tensor_p1_stiffness(
    double* A, const double* coordinate_dofs);
```

`A` is a contiguous, row-major `ndofs × ndofs` matrix. `coordinate_dofs`
contains the four vertex coordinates as a contiguous `4 × 3` array. The
kernel accumulates into `A`, so zero it before computing a new element matrix.
No memref descriptors, MLIR headers, Python, JIT or MLIR runtime are needed by
the executable.

## Build and run

Generate the kernel and its independent reference matrix from `demo/`:

```bash
(cd .. && python3 generate_kernel.py 1)
cmake -S . -B build_bare
cmake --build build_bare
./build_bare/run_kernel
```

If the LLVM tools are not on `PATH`, configure with their explicit paths:

```bash
cmake -S . -B build_bare \
  -DMLIR_OPT=$LLVM_BUILD_DIR/bin/mlir-opt \
  -DMLIR_TRANSLATE=$LLVM_BUILD_DIR/bin/mlir-translate \
  -DLLC=$LLVM_BUILD_DIR/bin/llc
```

For another degree, generate the corresponding files and use a separate build
folder, for example:

```bash
(cd .. && python3 generate_kernel.py 2)
cmake -S . -B build_bare_p2 -DKERNEL_DEGREE=2
cmake --build build_bare_p2
./build_bare_p2/run_kernel
```

CMake caches `KERNEL_MLIR`, `KERNEL_EXPECTED` and `MLIR_LOWERING_PIPELINE`.
Use a fresh build directory when switching degree or migrating an existing
build from the descriptor ABI. Alternatively, explicitly reset the relevant
cache entries. `KERNEL_DEGREE` selects both the C++ symbol name and matrix size;
the kernel and reference file must agree with it.

`run_kernel` repeats the call 10,000 times, clearing `A` before each call, then
checks the final matrix against the independently computed reference. A
successful check ends with `PASS`. This is a correctness check, not a timed
benchmark. P1, P2 and P3 have been compiled and run with the bare-pointer ABI.

## Compilation pipeline

```text
MLIR kernel -> mlir-opt -> mlir-translate -> llc -> kernel.o
C++ caller + kernel.o -> executable
```

CMake drives the lowering, translation, object compilation and linking. Kernel
generation remains a separate step because it needs the Python environment
with Basix and NumPy. The native build only needs the LLVM command-line tools
and a C++ compiler.

The lowering pipeline uses:

```text
convert-func-to-llvm{use-bare-ptr-memref-call-conv=true}
```

This turns each statically shaped memref argument into a single data pointer.
The shape and strides are known to the compiled kernel. The caller uses
`tabulate_tensor_p<degree>_stiffness` directly, without the `_mlir_ciface_`
prefix. The input's `llvm.emit_c_interface` attribute can remain: bare-pointer
lowering does not generate that descriptor wrapper.

## Other kernels

Bare-pointer lowering requires statically shaped memrefs with default layouts.
The generated stiffness kernels satisfy this for both output and coordinates.
The general UFLx emitter currently uses a dynamically sized coefficient memref;
coefficient-bearing kernels need a static per-cell coefficient-buffer size or
a separate wrapper before using this convention.

An alternative `KERNEL_MLIR` must match the caller's symbol, argument order,
shapes and reference data. Merely changing that path does not adapt the caller
to a different signature. Kernels containing math operations may also require
`convert-math-to-llvm` in `MLIR_LOWERING_PIPELINE`.

This change is local to the C++ demo's AOT pipeline. The Python harness retains
its existing descriptor calling convention.
