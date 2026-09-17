# Calling an AOT-compiled MLIR kernel from C++

A minimal, fully offline example of the ahead-of-time path discussed for
using a generated kernel from C++: save the kernel's MLIR text, lower and
compile it to a real object file entirely with command-line tools, then
link and call it from a plain C++ program with no MLIR/LLVM runtime, JIT,
or Python involved at run time.

This uses the coefficient-free P`<degree>` tetrahedron stiffness kernels
`generate_kernel.py` produces (`../kernels/p<degree>_stiffness.mlir`, for
whichever degree you generate) rather than one of the Coefficient-carrying
kernels from `generate_mlir_module`, specifically to keep the calling
convention down to two memref arguments instead of three. See `main.cpp`'s
docstring comment for exactly what changes to extend this to a
Coefficient-based (mass/stiffness "linear action") kernel -- it's one more
argument of the same shape, nothing structurally different.

Files: `CMakeLists.txt` (the build), `main.cpp` (the program that calls the kernel).

## Pipeline

```
form -> ... -> MLIR text            (generate_kernel.py / generate_mlir_module)
     -> mlir-opt                    (dialect lowering: memref/arith/scf/func -> llvm)
     -> mlir-translate              (llvm dialect -> LLVM IR)
     -> llc                         (LLVM IR -> native object file, kernel.o)
     -> g++/clang++                 (ordinary C++ compile + link against kernel.o)
```

Nothing past the first arrow needs Python, MLIR's Python bindings, or an
`ExecutionEngine` -- `mlir-opt`, `mlir-translate` and `llc` are the same
command-line tools your `llvm-project` build already produces (see the
package README's build instructions), and `kernel.o` is a completely
ordinary object file once `llc` is done with it.

## Build and run

Step 0, either way -- generate the kernel's MLIR text *and* its
independent reference matrix (needs only basix + numpy, no MLIR bindings,
so this stays a separate manual step rather than something CMake drives
itself -- see "Why CMake doesn't run generate_kernel.py" below):

```bash
(cd .. && python3 generate_kernel.py 1)
# writes ../kernels/p1_stiffness.mlir and ../kernels/p1_stiffness_expected.txt
```

Swap `1` for whatever degree you want to build against (`2`, `3`, ...) --
`generate_kernel.py` always writes both files together, under names
`main.cpp`/`CMakeLists.txt` expect to find as a matching pair.

### Build

```bash
cmake -B build
cmake --build build
./build/run_kernel
```

This builds and runs the P1 kernel by default. To build against a
different degree (once you've generated it in step 0):

```bash
cmake -B build_p2 -DKERNEL_DEGREE=2
cmake --build build_p2
./build_p2/run_kernel
```

`KERNEL_DEGREE` (default `1`) is a CMake *cache* variable, and it only
drives `KERNEL_MLIR`'s and `KERNEL_EXPECTED`'s *default* values, computed
once at first configure -- so switching degree needs either a fresh build
directory (as above) or `-DKERNEL_MLIR=...`/`-DKERNEL_EXPECTED=...`
passed explicitly alongside `-DKERNEL_DEGREE=...` in an existing one
(cache variables aren't recomputed from each other on reconfigure).
`KERNEL_DEGREE` is also passed straight to `main.cpp` as a compile
definition: it's what `main.cpp` token-pastes into the kernel's actual
linker symbol name (`tabulate_tensor_p<degree>_stiffness`) and uses to
compute `ndofs` at compile time, so it has to agree with whichever
`KERNEL_MLIR` you're actually linking against.

`CMakeLists.txt` drives the whole rest of the pipeline itself: an
`add_custom_command` runs `mlir-opt -> mlir-translate -> llc` to produce
`build/kernel.o` (re-running automatically if `KERNEL_MLIR` changes),
then links it straight into `run_kernel` alongside `main.cpp` --
`kernel.o` is registered as an `EXTERNAL_OBJECT`/`GENERATED` source, which
is the standard CMake idiom for feeding a build-time-generated object file
into a target as if it were any other translation unit.

If `mlir-opt`/`mlir-translate`/`llc` aren't on `PATH`, point CMake at your
`llvm-project` build directly:

```bash
cmake -B build \
  -DMLIR_OPT=$LLVM_BUILD_DIR/bin/mlir-opt \
  -DMLIR_TRANSLATE=$LLVM_BUILD_DIR/bin/mlir-translate \
  -DLLC=$LLVM_BUILD_DIR/bin/llc
```

`-DKERNEL_MLIR=/path/to/other.mlir` points the whole pipeline at a
different kernel entirely (eg one of `generate_mlir_module`'s own dumps
for a real UFLx form) instead of one of `generate_kernel.py`'s own
per-degree outputs -- if you do this, also override `-DKERNEL_EXPECTED=`
to a matching reference file (or drop `run_kernel`'s comparison from
`main.cpp` entirely if there isn't one), since a stale default
`KERNEL_EXPECTED` almost certainly won't have the same `ndofs` and
`main.cpp` will refuse to run rather than compare nonsense.

`build/run_kernel` is the compiled program; the three tool invocations
above are exactly what `CMakeLists.txt`'s `add_custom_command` runs, if
you'd rather see them run without CMake -- copy them out of
`CMakeLists.txt` directly (`build_kernel.sh` did this in an earlier version
of this demo; it's gone now, since keeping the exact same pipeline string
correct in two places was a bigger liability than the convenience was
worth once CMake did the same job).

Expected output ends with `PASS` -- `main.cpp` checks the kernel's result
against the `p<degree>_stiffness_expected.txt` reference matrix
`generate_kernel.write_expected_reference()` computed independently in
Python (read at run time, not pasted into the C++ file, so it stays in
sync automatically as you regenerate a kernel), so this doubles as a
correctness check of the whole AOT pipeline for whichever degree you
built, not just a demo that it runs. (Confirmed working end to end at P1
-- this is what prompted adding CMakeLists.txt in the first place.)

**What I verified here vs. couldn't:** this sandbox has `llvm-project`
checked out but not built, so `mlir-opt`/`mlir-translate`/`llc` aren't
available -- I couldn't run the real lowering or link a real `kernel.o`.
I did install `cmake` (via `pip install cmake`, no root available here)
and configure + build this `CMakeLists.txt` against three stub
`mlir-opt`/`mlir-translate`/`llc` scripts that just copy/touch fake output
files, purely to exercise the custom-command chain, the
`EXTERNAL_OBJECT`/`GENERATED` linking setup, and (for this
`KERNEL_DEGREE` update) the per-degree defaulting of `KERNEL_MLIR`/
`KERNEL_EXPECTED` and the `target_compile_definitions` plumbing into
`main.cpp`: for both `KERNEL_DEGREE=1` (the default) and `-DKERNEL_DEGREE=2`,
configure picked the right `p<degree>_stiffness.mlir`/`_expected.txt`
pair, the `kernel_obj` custom target ran all three tool invocations
against the right per-degree input file, `main.cpp` compiled with the
right `-DKERNEL_DEGREE=<N>` and `-DKERNEL_EXPECTED=...`, the link against
a *real* stand-in object file exposing the exact expected
`_mlir_ciface_tabulate_tensor_p<degree>_stiffness` symbol name succeeded
for both degrees, and reconfiguring with a degree that hasn't been
generated yet (`-DKERNEL_DEGREE=99`) failed configure with the clear
"generate it first" message rather than a confusing later error.
`main.cpp` itself (the token-pasting macros, `constexpr` ndofs, and
runtime `read_expected()`/mismatch-detection logic) was separately
compiled, linked and run standalone (outside CMake) for `KERNEL_DEGREE`
1, 2 and 3, including its two error paths (missing reference file; a
reference file whose `ndofs` doesn't match `KERNEL_DEGREE`) -- all
produced exactly the expected output. What's still untested is the real
`mlir-opt`/`mlir-translate`/`llc` toolchain's actual output, same caveat
as `main.cpp`'s `StridedMemRefType` layout before it.

## Why CMake doesn't run generate_kernel.py

`generate_kernel.py` needs `basix` (and, transitively through
`../harness.py`'s own imports if you look at more than just the functions
this demo uses, the MLIR Python bindings) -- a real Python environment
CMake has no reliable way to locate or activate on its own. Keeping that
step manual means `CMakeLists.txt` only has to depend on tools that are
trivial to `find_program()` (`mlir-opt`, `mlir-translate`, `llc`), and
`KERNEL_MLIR` stays a plain input file CMake's normal dependency tracking
(`DEPENDS` on the `add_custom_command`) already knows how to react to --
edit or regenerate the `.mlir` file and the next build re-lowers it
automatically.

## The calling convention (`_mlir_ciface_*`)

Every kernel `generate_mlir_module`/`generate_kernel.py` emits carries the
`llvm.emit_c_interface` attribute (visible right in the `.mlir` text, on
the `func.func` line). That attribute makes the lowering produce *two*
entry points: the kernel's "raw" LLVM-level function (memref arguments
unpacked into individual scalar/pointer fields -- this is what
`ExecutionEngine`'s packed-args calling convention in `../harness.py`
uses, not what you want from plain C/C++), and a friendlier
`_mlir_ciface_<kernel_name>` wrapper built specifically for C callers,
which takes one pointer per memref argument, each pointing at a small
descriptor struct: base pointer, data pointer, offset, `sizes[rank]`,
`strides[rank]`.

`main.cpp` reproduces that struct by hand (`StridedMemRefType<T, Rank>`)
so the example needs no MLIR include path to build. The authoritative
version is `mlir::StridedMemRefType<T, Rank>` in
`mlir/ExecutionEngine/CRunnerUtils.h` in your LLVM checkout -- a small,
dependency-free header -- and `#include`-ing it directly instead is a
perfectly good alternative to hand-duplicating the struct if you'd rather
not maintain that copy yourself.

`main.cpp` also builds `KERNEL_CIFACE_NAME` (the exact symbol it declares
`extern "C"` and calls) from `KERNEL_DEGREE` via preprocessor
token-pasting, since the kernel's actual name -- and hence its
`_mlir_ciface_*` wrapper's name -- includes the degree
(`tabulate_tensor_p<degree>_stiffness`, matching
`generate_kernel.py`'s own naming exactly). If `KERNEL_DEGREE` doesn't
match the degree `KERNEL_MLIR` was generated for, the link step fails
with an undefined-symbol error naming the mismatched symbol, rather than
silently calling the wrong kernel.

## Why no runtime library

`../harness.py` builds its `ExecutionEngine` with `shared_libs=[]` for
every CPU kernel in this repo, with a comment noting that
`libmlir_runner_utils`/`libmlir_c_runner_utils` would only be needed if a
kernel called into printing or memref-allocation helpers -- none of the
stiffness/mass kernels here do. That's why `g++ ... main.cpp kernel.o` is
a complete link line: `kernel.o` is self-contained native code, nothing
else to add beyond libc/libm.
