"""
JIT harness for the P1 stiffness kernel.

Prereqs (after building LLVM/MLIR with Python bindings enabled -- see the
package README's "Building LLVM/MLIR with Python bindings" section):
    export PYTHONPATH=/path/to/llvm-project/build/tools/mlir/python_packages/mlir_core:$PYTHONPATH
    python3 -c "import mlir; print(mlir.__file__)"   # sanity check

Run (from the uflx_mlir package root):
    python3 demo/harness.py
"""

import ctypes
import time
from pathlib import Path

import numpy as np

from mlir.ir import Context, Module
from mlir.passmanager import PassManager
from mlir.execution_engine import ExecutionEngine
from mlir.runtime import get_ranked_memref_descriptor

KERNEL_PATH = Path(__file__).parent / "kernels" / "p1_stiffness.mlir"
KERNEL_NAME = "tabulate_tensor_p1_stiffness"

# Lowering pipeline: memref/arith/func -> llvm dialect.
# Pass names occasionally change between LLVM releases -- if PassManager.parse
# raises "unknown pass", check `mlir-opt --help` from your build for the
# current name of the memref-to-llvm and func-to-llvm passes.
PIPELINE = (
    "builtin.module("
    "convert-arith-to-llvm,"
    "finalize-memref-to-llvm,"
    "convert-func-to-llvm,"
    "reconcile-unrealized-casts"
    ")"
)

# Pipeline for kernels that use scf.for (all quadrature-loop kernels from
# generate_kernel.py, including P1) -- structured control flow needs lowering to cf
# before it can reach llvm. Pass names confirmed via `mlir-opt --help` on
# this build (LLVM release/18.x), not guessed.
QUADRATURE_PIPELINE = (
    "builtin.module("
    "convert-scf-to-cf,"
    "convert-arith-to-llvm,"
    "finalize-memref-to-llvm,"
    "convert-cf-to-llvm,"
    "convert-func-to-llvm,"
    "reconcile-unrealized-casts"
    ")"
)

# Same as QUADRATURE_PIPELINE, plus convert-math-to-llvm -- needed by any
# kernel using ops from the `math` dialect (e.g. math.absf, emitted for
# abs(detJ) by the UFLx-driven generator in ../uflx_mlir/emit.py, since
# uflx.geometry.JacobianDeterminant always takes abs() to handle either
# cell orientation). Harmless to run on a kernel with no math.* ops in it,
# so this could have just been folded into QUADRATURE_PIPELINE above, but
# keeping it separate avoids touching a pipeline the P2-P5 kernels already
# validate against.
UFLX_PIPELINE = (
    "builtin.module("
    "convert-scf-to-cf,"
    "convert-math-to-llvm,"
    "convert-arith-to-llvm,"
    "finalize-memref-to-llvm,"
    "convert-cf-to-llvm,"
    "convert-func-to-llvm,"
    "reconcile-unrealized-casts"
    ")"
)

# Same as UFLX_PIPELINE, but with canonicalize/cse/loop-invariant-code-motion
# run FIRST, while loops are still structured `scf.for` -- added after
# ffcx_compare_uflx.py showed the UFLx-generated P3 stiffness kernel running
# ~27x slower per call than FFCx's. The generated loop nest is i(20) / j(20)
# / quadrature-point(14), with the WHOLE per-entry expression -- including
# the Jacobian/detJ/cofactor geometry terms, which only depend on the cell's
# coordinates and are invariant across all 5600 (i,j,q) combinations --
# flattened into the innermost loop body (confirmed by inspecting a
# generated p3_stiffness_uflx.mlir: a single ~1400-op basic block with no
# per-node reuse across loop levels). generate_kernel.py's hand-written
# kernel avoids this entirely by looping quadrature outermost and caching
# each dof's physical gradient in a scratch array once per quadrature point.
# ../uflx_mlir/emit.py's generate_mlir_module now does this same
# loop-hoisting/fission restructuring itself (see ../uflx_mlir/hoist.py),
# which is why UFLX_OPTIMIZED_PIPELINE below is no longer needed to close
# the gap -- kept here for historical comparison (the "baseline" row in
# ffcx_compare_uflx.py's output) and because it's still a useful sanity
# check that MLIR's own LICM doesn't already do hoist.py's job for free.
UFLX_OPTIMIZED_PIPELINE = (
    "builtin.module("
    "canonicalize,"
    "cse,"
    "loop-invariant-code-motion,"
    "convert-scf-to-cf,"
    "convert-math-to-llvm,"
    "convert-arith-to-llvm,"
    "finalize-memref-to-llvm,"
    "convert-cf-to-llvm,"
    "convert-func-to-llvm,"
    "reconcile-unrealized-casts"
    ")"
)

# Cache of ExecutionEngine.lookup() results, keyed by id(engine) -- see
# build_caller() below. A plain dict is fine for a prototype with a handful of
# long-lived engines; it would leak if engines were churned rapidly and their
# ids got reused, which doesn't happen here.
_caller_cache = {}


def build_engine_from_text(
    mlir_text: str, pipeline: str = PIPELINE, opt_level: int = 3
) -> ExecutionEngine:
    """Core primitive behind build_engine_from() below -- takes MLIR text
    directly instead of a file path, e.g. for timing a full
    form -> text -> engine pipeline without needing a temp file (see
    ffcx_compare_uflx.py, which needs generate_mlir_module()'s
    str(module)-serialized text JIT'd without writing it to kernels/*.mlir
    first).

    opt_level is passed straight to ExecutionEngine, which forwards it to
    LLVM's own optimizer on the lowered IR (separate from any MLIR-level
    passes in `pipeline`) -- e.g. loop-invariant code motion, vectorization.
    Explicit 3 here rather than relying on whatever ExecutionEngine's default
    is, since FFCx's cffi-compiled C code gets real optimization by default
    via the system compiler and this should get a fair comparison too."""
    with Context():
        module = Module.parse(mlir_text)
        pm = PassManager.parse(pipeline)
        pm.run(module.operation)
        # shared_libs=[] here; add libmlir_runner_utils/libmlir_c_runner_utils
        # if you later use printing/memref helper calls from the kernel itself.
        engine = ExecutionEngine(module, opt_level=opt_level)
        return engine


def build_engine_from(
    kernel_path: Path, pipeline: str = PIPELINE, opt_level: int = 3
) -> ExecutionEngine:
    """General entry point: build an ExecutionEngine for any kernel file /
    pipeline pair. build_engine() below is the P1-specific convenience
    wrapper kept for backward compatibility. Thin wrapper around
    build_engine_from_text() -- see that docstring for the opt_level note."""
    return build_engine_from_text(kernel_path.read_text(), pipeline, opt_level)


def build_engine() -> ExecutionEngine:
    return build_engine_from(KERNEL_PATH, QUADRATURE_PIPELINE)


def build_engine_from_module(
    module: Module, pipeline: str = PIPELINE, opt_level: int = 3
) -> ExecutionEngine:
    """Like build_engine_from, but for a Module already built in memory via
    the Python op-builder API (see ../uflx_mlir/emit.py) instead of parsed
    from MLIR text -- skips the str(module) -> Module.parse() round-trip
    entirely, running the same pass pipeline and JIT directly on the module
    the caller already built (and, in uflx_mlir.emit.generate_mlir_module's
    case, already verified with module.operation.verify()).

    module.context is the Context it was created under -- Module.create()
    keeps a reference alive, so this doesn't need a fresh Context() the way
    build_engine_from() does when parsing text from scratch. PassManager and
    ExecutionEngine both need that same context current, hence the `with`.
    """
    with module.context:
        pm = PassManager.parse(pipeline)
        pm.run(module.operation)
        engine = ExecutionEngine(module, opt_level=opt_level)
        return engine


def check_lowering(kernel_path: Path, pipeline: str = PIPELINE) -> str:
    """Parse a kernel and run the given pass pipeline WITHOUT building an
    ExecutionEngine or JIT'ing anything -- the same "does this lower cleanly"
    check `mlir-opt <file> --pass-pipeline=...` does from the shell, but via
    the same Python bindings the rest of this file uses, so there's no
    separate binary invocation needed. A bad pipeline or invalid IR raises a
    normal Python exception here; this step was never the crash risk (that
    was ExecutionEngine.lookup()'s calling convention, a different step).
    Returns the lowered module's text."""
    mlir_text = kernel_path.read_text()
    with Context():
        module = Module.parse(mlir_text)
        pm = PassManager.parse(pipeline)
        pm.run(module.operation)
        return str(module)


def build_caller(engine: ExecutionEngine):
    """Look up the kernel's JIT'd entry point once and return a closure that
    calls it directly, bypassing ExecutionEngine.invoke()'s per-call symbol
    lookup and argument-list bookkeeping -- measured ~15x faster per call for
    this kernel (see ffcx_compare.py, which found this the hard way: three
    segfaults before confirming the real calling convention).

    engine.lookup(name) returns a bound ctypes callable using MLIR's generic
    "packed" convention: ONE argument, a pointer to an array of void*, where
    slot i holds the ADDRESS of the i-th real argument (confirmed via
    introspection: its .argtypes is (ctypes.c_void_p,) -- a single arg, not
    one-per-kernel-argument). This is exactly what invoke() must assemble
    internally from the args you pass it.
    """
    return build_caller_for(engine, KERNEL_NAME)


def build_caller_for(engine: ExecutionEngine, kernel_name: str):
    """General entry point behind build_caller() -- see its docstring."""
    raw_fn = engine.lookup(kernel_name)

    def call(A: np.ndarray, coords: np.ndarray) -> None:
        A_desc_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(A)))
        coords_desc_pp = ctypes.pointer(
            ctypes.pointer(get_ranked_memref_descriptor(coords))
        )
        packed = (ctypes.c_void_p * 2)(
            ctypes.cast(A_desc_pp, ctypes.c_void_p).value,
            ctypes.cast(coords_desc_pp, ctypes.c_void_p).value,
        )
        raw_fn(packed)

    return call


def run_kernel(engine: ExecutionEngine, coords: np.ndarray) -> np.ndarray:
    assert coords.shape == (4, 3) and coords.dtype == np.float64
    A = np.zeros((4, 4), dtype=np.float64)

    call = _caller_cache.get(id(engine))
    if call is None:
        call = build_caller(engine)
        _caller_cache[id(engine)] = call
    call(A, coords)
    return A


def reference_p1_stiffness(coords: np.ndarray) -> np.ndarray:
    """Direct numpy computation for validation -- same math, independent path.

    coords: 4x3 array of tetrahedron vertex coordinates (x0..x3)."""
    x0, x1, x2, x3 = coords
    J = np.column_stack([x1 - x0, x2 - x0, x3 - x0])  # 3x3 Jacobian
    detJ = np.linalg.det(J)
    Jinv = np.linalg.inv(J)
    ref_grads = np.array(
        [[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    phys_grads = ref_grads @ Jinv  # each row: physical grad of phi_i
    volume = detJ / 6.0
    return (phys_grads @ phys_grads.T) * volume


def main():
    coords = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )

    engine = build_engine()
    A = run_kernel(engine, coords)
    A_ref = reference_p1_stiffness(coords)

    print("MLIR JIT result:\n", A)
    print("Reference result:\n", A_ref)
    np.testing.assert_allclose(A, A_ref, rtol=1e-12, atol=1e-12)
    print("MATCH: within 1e-12")

    # Rough kernel-call timing (compile time is a separate, one-off cost above).
    n = 200_000
    t0 = time.perf_counter()
    for _ in range(n):
        run_kernel(engine, coords)
    t1 = time.perf_counter()
    print(f"{n} invocations in {t1 - t0:.4f}s -> {(t1 - t0) / n * 1e6:.3f} us/call")
    print("Note: run_kernel() still builds a fresh memref descriptor per call")
    print("(coords/A can differ between calls); only the engine.lookup() -> raw")
    print("function pointer step is cached. See ffcx_compare.py for a stricter")
    print("apples-to-apples timing where descriptors are also built once and reused.")


if __name__ == "__main__":
    main()
