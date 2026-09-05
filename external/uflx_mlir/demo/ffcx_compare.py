"""
Side-by-side FFCx-vs-MLIR comparison for a CG-`degree` Laplacian stiffness
kernel: every degree uses a basix-generated quadrature-loop kernel from
generate_kernel.py (run that first).

Measures:
  - cold "compile" time (FFCx: codegen + native compile; MLIR: parse + lower +
    ExecutionEngine build)
  - steady-state per-call time, using a minimal-overhead wrapper on both sides
    (buffers/descriptors built once, reused across all calls -- isolates
    invoke+kernel cost from Python-side allocation/marshaling noise)
  - a combined "compile once, then call N times" cold-start total

No dolfinx/PETSc/MPI dependency -- only ufl + basix + ffcx for the FFCx side.

Run (same venv, fenics-ffcx installed; generate_kernel.py <degree> first):
    python3 demo/ffcx_compare.py        # degree 1 (default)
    python3 demo/ffcx_compare.py 2
    python3 demo/ffcx_compare.py 3
"""

import ctypes
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

import basix.ufl
import ufl
from ffcx.codegeneration.jit import compile_forms
from mlir.runtime import get_ranked_memref_descriptor

sys.path.insert(0, str(Path(__file__).parent))
import harness as mlir_harness

N_CALLS_DEFAULT = 5000  # steady-state repetitions -- large enough that per-call
# noise (scheduling jitter, etc.) averages out; the mean is what gets reported,
# and its variance shrinks like 1/n_calls, so more samples = less noise in the
# printed us/call number. Override with a second CLI arg if this is too slow
# (FFCx's compile step happens once regardless, so N_CALLS only affects the
# steady-state loop's wall-clock, not the cold-compile numbers).

# Scalene tetrahedron -- avoids masking an index/transpose bug behind
# accidental symmetry (same coordinates run_higher_order.py and
# generate_kernel.py's own cross-check use).
COORDS = np.array(
    [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
    dtype=np.float64,
)


def ndofs_for_degree(degree: int) -> int:
    return (degree + 1) * (degree + 2) * (degree + 3) // 6  # tet Lagrange DOF count


def build_form(degree: int):
    # Must match generate_kernel.py's LAGRANGE_VARIANT exactly -- for degree<=2
    # there's only one point per edge (or none), so any variant coincides and
    # this wouldn't have shown up; degree 3 has two points per edge, where
    # equispaced (1/3, 2/3 along the edge) and DOLFINx's usual default variant
    # genuinely differ, which is exactly what made FFCx and MLIR disagree here.
    element = basix.ufl.element(
        "Lagrange", "tetrahedron", degree,
        lagrange_variant=basix.LagrangeVariant.equispaced,
    )
    coord_element = basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,))
    domain = ufl.Mesh(coord_element)
    V = ufl.FunctionSpace(domain, element)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    return ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx


def ffcx_compile(degree: int):
    """Returns (compile_seconds, fast_call, A). fast_call() reruns the kernel
    in place into A via a minimal-overhead cffi wrapper (buffers built once)."""
    a = build_form(degree)
    ndofs = ndofs_for_degree(degree)
    cache_dir = Path(tempfile.mkdtemp(prefix="ffcx_jit_"))
    t0 = time.perf_counter()
    ufcx_forms, module, code = compile_forms([a], options={}, cache_dir=cache_dir)
    t1 = time.perf_counter()

    ffi = module.ffi
    ufcx_form = ufcx_forms[0]
    offsets = [ufcx_form.form_integral_offsets[i] for i in range(4)]
    start, end = offsets[0], offsets[1]
    if end - start < 1:
        raise RuntimeError(f"No cell integrals found (offsets={offsets})")
    tabulate_tensor = ufcx_form.form_integrals[start].tabulate_tensor_float64

    # gdim is already 3 for a tetrahedron -- no stride-3 z-padding needed
    # (that was only for the 2D triangle case, where ufcx's coordinate_dofs
    # convention still pads to 3 components per vertex).
    coords_flat = np.ascontiguousarray(COORDS.reshape(-1))
    entity_local_index = np.zeros(1, dtype=np.int32)
    quadrature_permutation = np.zeros(1, dtype=np.uint8)
    A = np.zeros((ndofs, ndofs), dtype=np.float64)

    A_buf = ffi.from_buffer("double[]", A)
    coords_buf = ffi.from_buffer("double[]", coords_flat)
    eli_buf = ffi.from_buffer("int[]", entity_local_index)
    qp_buf = ffi.from_buffer("uint8_t[]", quadrature_permutation)

    def fast_call():
        tabulate_tensor(A_buf, ffi.NULL, ffi.NULL, coords_buf, eli_buf, qp_buf, ffi.NULL)

    shutil.rmtree(cache_dir, ignore_errors=True)
    return t1 - t0, fast_call, A


def mlir_compile(degree: int):
    """Returns (compile_seconds, fast_calls, A). fast_calls is a dict of two
    variants, both reusing pre-built descriptors:
      'invoke'        -- via ExecutionEngine.invoke()
      'direct_ctypes' -- via ExecutionEngine.lookup() + the packed-args
                          calling convention (bypasses invoke()'s per-call
                          overhead; see harness.build_caller_for)."""
    kernel_path = Path(__file__).parent / "kernels" / f"p{degree}_stiffness.mlir"
    if not kernel_path.exists():
        hint = f" -- run `python3 demo/generate_kernel.py {degree}` first"
        raise FileNotFoundError(f"{kernel_path} doesn't exist{hint}")
    kernel_name = f"tabulate_tensor_p{degree}_stiffness"
    pipeline = mlir_harness.QUADRATURE_PIPELINE
    ndofs = ndofs_for_degree(degree)

    t0 = time.perf_counter()
    engine = mlir_harness.build_engine_from(kernel_path, pipeline)
    t1 = time.perf_counter()

    A = np.zeros((ndofs, ndofs), dtype=np.float64)
    A_desc = get_ranked_memref_descriptor(A)
    coords_desc = get_ranked_memref_descriptor(COORDS)

    A_desc_pp = ctypes.pointer(ctypes.pointer(A_desc))
    coords_desc_pp = ctypes.pointer(ctypes.pointer(coords_desc))

    def call_invoke():
        engine.invoke(kernel_name, A_desc_pp, coords_desc_pp)

    fast_calls = {"invoke": call_invoke}

    try:
        addr = engine.lookup(kernel_name)
        packed = (ctypes.c_void_p * 2)(
            ctypes.cast(A_desc_pp, ctypes.c_void_p).value,
            ctypes.cast(coords_desc_pp, ctypes.c_void_p).value,
        )

        def call_direct():
            addr(packed)

        fast_calls["direct_ctypes"] = call_direct
    except Exception as e:
        print(f"(direct ctypes lookup path unavailable: {e!r} -- skipping that variant)", flush=True)

    return t1 - t0, fast_calls, A


def time_calls(fast_call, n):
    t0 = time.perf_counter()
    for _ in range(n):
        fast_call()
    t1 = time.perf_counter()
    return t1 - t0


def reference_for(degree: int, coords: np.ndarray) -> np.ndarray:
    if degree == 1:
        return mlir_harness.reference_p1_stiffness(coords)
    import generate_kernel  # only needed for degree > 1

    return generate_kernel.reference_stiffness(coords, degree)


def main():
    degree = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n_calls = int(sys.argv[2]) if len(sys.argv) > 2 else N_CALLS_DEFAULT
    print(f"=== degree {degree}, N_CALLS={n_calls} ===", flush=True)
    A_ref = reference_for(degree, COORDS)

    print("--- FFCx ---", flush=True)
    ffcx_compile_s, ffcx_call, A_ffcx = ffcx_compile(degree)
    ffcx_call()
    np.testing.assert_allclose(A_ffcx, A_ref, rtol=1e-9, atol=1e-8)
    print(f"compile: {ffcx_compile_s * 1e3:.3f} ms  (validated: MATCH)", flush=True)
    ffcx_calls_s = time_calls(ffcx_call, n_calls)
    ffcx_us_per_call = ffcx_calls_s / n_calls * 1e6
    print(f"{n_calls} calls: {ffcx_calls_s * 1e3:.3f} ms -> {ffcx_us_per_call:.3f} us/call", flush=True)
    ffcx_total_s = ffcx_compile_s + ffcx_calls_s
    print(f"compile + {n_calls} calls, cold start: {ffcx_total_s * 1e3:.3f} ms", flush=True)

    print("\n--- MLIR ---", flush=True)
    mlir_compile_s, mlir_calls, A_mlir = mlir_compile(degree)
    print(f"compile: {mlir_compile_s * 1e3:.3f} ms", flush=True)

    mlir_results = {}
    for variant, fast_call in mlir_calls.items():
        fast_call()
        np.testing.assert_allclose(A_mlir, A_ref, rtol=1e-9, atol=1e-8)
        calls_s = time_calls(fast_call, n_calls)
        us_per_call = calls_s / n_calls * 1e6
        total_s = mlir_compile_s + calls_s
        print(f"[{variant}] validated: MATCH -- {n_calls} calls: {calls_s * 1e3:.3f} ms "
              f"-> {us_per_call:.3f} us/call; compile + calls cold start: {total_s * 1e3:.3f} ms",
              flush=True)
        mlir_results[variant] = (us_per_call, total_s)

    best_variant = min(mlir_results, key=lambda k: mlir_results[k][0])
    mlir_us_per_call, mlir_total_s = mlir_results[best_variant]

    print(f"\n--- Summary (degree={degree}, N={n_calls}, MLIR best variant: {best_variant}) ---",
          flush=True)
    print(f"compile time,  FFCx / MLIR: {ffcx_compile_s / mlir_compile_s:.3f}x", flush=True)
    print(f"per-call time, FFCx / MLIR: {ffcx_us_per_call / mlir_us_per_call:.3f}x", flush=True)
    print(f"cold total,    FFCx / MLIR: {ffcx_total_s / mlir_total_s:.3f}x", flush=True)
    print(flush=True)
    print(f"Both call paths reuse pre-built buffers/descriptors across all {n_calls}", flush=True)
    print("calls (built once, not rebuilt per call), so per-call numbers isolate", flush=True)
    print("invoke+kernel cost rather than Python-side allocation/marshaling noise.", flush=True)


if __name__ == "__main__":
    main()
