"""Full-pipeline comparison, for a Coefficient's linear action: UFLx form ->
MLIR -> JIT'd kernel, versus legacy-UFL form -> FFCx -> compiled kernel, for
`inner(grad(w), grad(v)) * dx` where `w` is a Coefficient (not a
TrialFunction) on a CG-`degree` tetrahedron space.

This is ffcx_compare_uflx.py's sibling for the Coefficient/linear-action
case rather than the bilinear (matrix-assembly) case: `w` is a genuine
Coefficient, so the local tensor is a length-ndofs VECTOR (the stiffness
matrix's action on w's own dof vector), not an ndofs x ndofs matrix, and
both sides' generated kernels take a third "coefficients" buffer argument
alongside the local tensor and the cell coordinates. See
uflx_codegeneration.algorithms.coefficients.insert_coefficient_functions
and uflx_mlir.emit's FunctionCall handling / _emit_coefficient_function for
how the UFLx -> MLIR side generates that: each distinct coefficient
(derivative, component) combination becomes its own small func.func doing
an scf.for dof-summation, called from the main kernel wherever the
coefficient's value is needed.

Same "compile" and steady-state per-call timing methodology as
ffcx_compare_uflx.py (see that file's docstring): compile_seconds covers
the WHOLE UFLx side (form -> quadrature/geometry/tabulation lowering ->
MLIR generation -> JIT) for a fair comparison against FFCx's
compile_forms() (UFL analysis + C codegen + a real C compile), and
per-call numbers reuse pre-built buffers/descriptors on both sides
throughout.

Run (same venv, fenics-ffcx installed):
    python3 demo/ffcx_compare_linear_action.py        # degree 1 (default)
    python3 demo/ffcx_compare_linear_action.py 2
    python3 demo/ffcx_compare_linear_action.py 3 50000  # override N_CALLS
"""

import ctypes
import shutil
import sys
import tempfile
import time
from pathlib import Path

import basix
import basix.ufl
import numpy as np
import ufl
from ffcx.codegeneration.jit import compile_forms
from mlir.runtime import get_ranked_memref_descriptor

sys.path.insert(0, str(Path(__file__).parent))
import generate_kernel
import harness as mlir_harness
from basix_uflx import element
from uflx import Coefficient, TestFunction, coordinate_element, dx, function_space, grad, inner

from uflx_mlir.emit import generate_mlir_module

N_CALLS_DEFAULT = 20_000  # see ffcx_compare.py's docstring for the reasoning

# Scalene tetrahedron -- same coordinates every other script in this repo
# uses, so results are directly comparable across files.
COORDS = np.array(
    [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
    dtype=np.float64,
)


def ndofs_for_degree(degree: int) -> int:
    return (degree + 1) * (degree + 2) * (degree + 3) // 6  # tet Lagrange DOF count


def build_coefficient_gradient_form(degree: int):
    """`inner(grad(w), grad(v)) * dx` for a Coefficient w on a P{degree}
    Lagrange tetrahedron space (coordinate map always P1/affine, same as
    ffcx_compare_uflx.py's build_stiffness_form). Returns (form, ndofs).
    """
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    w = Coefficient(space)
    v = TestFunction(space)
    return inner(grad(w), grad(v)) * dx, e.dim


def build_legacy_ufl_form(degree: int):
    """The same form as build_coefficient_gradient_form, but expressed in
    legacy UFL -- what FFCx actually consumes. Must use the same Lagrange
    variant (equispaced) as the UFLx side -- see ffcx_compare_uflx.py's
    build_legacy_ufl_form docstring for why this matters starting at
    degree 3.
    """
    ufl_element = basix.ufl.element(
        "Lagrange",
        "tetrahedron",
        degree,
        lagrange_variant=basix.LagrangeVariant.equispaced,
    )
    coord_element = basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,))
    domain = ufl.Mesh(coord_element)
    V = ufl.FunctionSpace(domain, ufl_element)
    w = ufl.Coefficient(V)
    v = ufl.TestFunction(V)
    return ufl.inner(ufl.grad(w), ufl.grad(v)) * ufl.dx


def ffcx_compile(degree: int, w_dofs: np.ndarray):
    """Returns (compile_seconds, fast_call, vec). Same shape as
    ffcx_compare_uflx.py's ffcx_compile(), except the local tensor is a
    length-ndofs vector (one Argument, not two) and a real coefficients
    buffer is passed instead of ffi.NULL -- FFCx's tabulate_tensor argument
    order is (A, w, constants, coordinate_dofs, entity_local_index,
    quadrature_permutation, custom_data); see ffcx_compare.py for the same
    convention with w left NULL.
    """
    a = build_legacy_ufl_form(degree)
    ndofs = ndofs_for_degree(degree)
    cache_dir = Path(tempfile.mkdtemp(prefix="ffcx_jit_"))
    t0 = time.perf_counter()
    ufcx_forms, module, _code = compile_forms([a], options={}, cache_dir=cache_dir)
    t1 = time.perf_counter()

    ffi = module.ffi
    ufcx_form = ufcx_forms[0]
    offsets = [ufcx_form.form_integral_offsets[i] for i in range(4)]
    start, end = offsets[0], offsets[1]
    if end - start < 1:
        raise RuntimeError(f"No cell integrals found (offsets={offsets})")
    tabulate_tensor = ufcx_form.form_integrals[start].tabulate_tensor_float64

    coords_flat = np.ascontiguousarray(COORDS.reshape(-1))
    entity_local_index = np.zeros(1, dtype=np.int32)
    quadrature_permutation = np.zeros(1, dtype=np.uint8)
    vec = np.zeros(ndofs, dtype=np.float64)
    w = np.ascontiguousarray(w_dofs, dtype=np.float64)

    vec_buf = ffi.from_buffer("double[]", vec)
    w_buf = ffi.from_buffer("double[]", w)
    coords_buf = ffi.from_buffer("double[]", coords_flat)
    eli_buf = ffi.from_buffer("int[]", entity_local_index)
    qp_buf = ffi.from_buffer("uint8_t[]", quadrature_permutation)

    def fast_call():
        tabulate_tensor(vec_buf, w_buf, ffi.NULL, coords_buf, eli_buf, qp_buf, ffi.NULL)

    shutil.rmtree(cache_dir, ignore_errors=True)
    return t1 - t0, fast_call, vec


def _mlir_fast_calls(engine, kernel_name: str, ndofs: int, w_dofs: np.ndarray):
    """Like ffcx_compare_uflx.py's _mlir_fast_calls (same dual 'invoke' /
    'direct_ctypes' variant pattern), but for a kernel whose form uses a
    Coefficient: generate_mlir_module only gives the assembly function a
    third memref<?xf64> argument in that case (see ../uflx_mlir/emit.py's
    own docstring and _OpCtx.coeffs_val), so the packed calling convention
    here has three descriptors instead of two.
    """
    vec = np.zeros(ndofs, dtype=np.float64)
    w = np.ascontiguousarray(w_dofs, dtype=np.float64)
    vec_desc = get_ranked_memref_descriptor(vec)
    coords_desc = get_ranked_memref_descriptor(COORDS)
    w_desc = get_ranked_memref_descriptor(w)

    vec_desc_pp = ctypes.pointer(ctypes.pointer(vec_desc))
    coords_desc_pp = ctypes.pointer(ctypes.pointer(coords_desc))
    w_desc_pp = ctypes.pointer(ctypes.pointer(w_desc))

    def call_invoke():
        engine.invoke(kernel_name, vec_desc_pp, coords_desc_pp, w_desc_pp)

    fast_calls = {"invoke": call_invoke}

    try:
        addr = engine.lookup(kernel_name)
        packed = (ctypes.c_void_p * 3)(
            ctypes.cast(vec_desc_pp, ctypes.c_void_p).value,
            ctypes.cast(coords_desc_pp, ctypes.c_void_p).value,
            ctypes.cast(w_desc_pp, ctypes.c_void_p).value,
        )

        def call_direct():
            addr(packed)

        fast_calls["direct_ctypes"] = call_direct
    except Exception as e:
        print(
            f"(direct ctypes lookup path unavailable: {e!r} -- skipping that variant)", flush=True
        )

    return fast_calls, vec


def uflx_compile(degree: int, w_dofs: np.ndarray, pipeline: str | None = None):
    """Full pipeline, direct-from-module: UFLx form -> generate_mlir_module's
    in-memory Module -> JIT, with no MLIR-text round-trip at all -- see
    ffcx_compare_uflx.py's uflx_compile docstring for the "compile" timing
    convention this mirrors exactly. pipeline defaults to
    mlir_harness.UFLX_PIPELINE; pass mlir_harness.UFLX_OPTIMIZED_PIPELINE to
    additionally run canonicalize/cse/loop-invariant-code-motion first.
    """
    pipeline = pipeline if pipeline is not None else mlir_harness.UFLX_PIPELINE
    form, ndofs = build_coefficient_gradient_form(degree)
    kernel_name = f"tabulate_tensor_p{degree}_coefficient_gradient_uflx"

    t0 = time.perf_counter()
    module = generate_mlir_module(
        form,
        degree=degree,
        kernel_name=kernel_name,
        cell=basix.CellType.tetrahedron,
        inline_geometry=True,
    )
    engine = mlir_harness.build_engine_from_module(module, pipeline)
    t1 = time.perf_counter()

    fast_calls, vec = _mlir_fast_calls(engine, kernel_name, ndofs, w_dofs)
    return t1 - t0, fast_calls, vec


def time_calls(fast_call, n):
    t0 = time.perf_counter()
    for _ in range(n):
        fast_call()
    t1 = time.perf_counter()
    return t1 - t0


def run_variant(label: str, compile_fn, degree: int, n_calls: int, vec_ref: np.ndarray):
    """Runs one 'compile, validate, time steady-state calls' pass -- shared
    by FFCx and the UFLx compile path below. Returns
    (compile_seconds, best_us_per_call, best_total_seconds, best_variant_name).
    """
    print(f"--- {label} ---", flush=True)
    compile_s, fast_calls, vec = compile_fn(degree)
    print(f"compile: {compile_s * 1e3:.3f} ms", flush=True)

    results = {}
    for variant, fast_call in fast_calls.items():
        # vec is shared across every variant in fast_calls (see
        # _mlir_fast_calls) -- the kernel is pure accumulate (matching
        # FFCx's/UFC's tabulate_tensor convention), so a variant running
        # after another has already accumulated into vec would otherwise
        # fail this check against a stale, non-zero buffer. Reset here,
        # once per variant, entirely outside the timed region below --
        # time_calls itself still reuses vec across all n_calls without
        # rezeroing, same as ffcx_compare_uflx.py.
        vec.fill(0.0)
        fast_call()
        np.testing.assert_allclose(vec, vec_ref, rtol=1e-9, atol=1e-8)
        calls_s = time_calls(fast_call, n_calls)
        us_per_call = calls_s / n_calls * 1e6
        total_s = compile_s + calls_s
        print(
            f"[{variant}] validated: MATCH -- {n_calls} calls: {calls_s * 1e3:.3f} ms "
            f"-> {us_per_call:.3f} us/call; compile + calls cold start: {total_s * 1e3:.3f} ms",
            flush=True,
        )
        results[variant] = (us_per_call, total_s)

    best = min(results, key=lambda k: results[k][0])
    return compile_s, results[best][0], results[best][1], best


# Same baseline/licm pipeline comparison as ffcx_compare_uflx.py -- see that
# file's PIPELINE_VARIANTS docstring for why UFLX_OPTIMIZED_PIPELINE exists
# (kept for historical comparison; ../uflx_mlir/hoist.py's own fission/
# hoisting analysis is what actually closes the gap now, MLIR's generic LICM
# doesn't rediscover it).
PIPELINE_VARIANTS = [
    ("baseline", None),  # None -> uflx_compile's own UFLX_PIPELINE default
    ("licm", "UFLX_OPTIMIZED_PIPELINE"),  # resolved via getattr below
]


def main():
    degree = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n_calls = int(sys.argv[2]) if len(sys.argv) > 2 else N_CALLS_DEFAULT
    print(f"=== degree {degree}, N_CALLS={n_calls} ===", flush=True)

    ndofs = ndofs_for_degree(degree)
    rng = np.random.default_rng(0)
    w_dofs = rng.standard_normal(ndofs)

    stiffness_ref = generate_kernel.reference_stiffness(COORDS, degree)
    vec_ref = stiffness_ref @ w_dofs

    print("--- FFCx (legacy UFL) ---", flush=True)
    ffcx_compile_s, ffcx_call, vec_ffcx = ffcx_compile(degree, w_dofs)
    ffcx_call()
    np.testing.assert_allclose(vec_ffcx, vec_ref, rtol=1e-9, atol=1e-8)
    print(f"compile: {ffcx_compile_s * 1e3:.3f} ms  (validated: MATCH)", flush=True)
    ffcx_calls_s = time_calls(ffcx_call, n_calls)
    ffcx_us_per_call = ffcx_calls_s / n_calls * 1e6
    print(
        f"{n_calls} calls: {ffcx_calls_s * 1e3:.3f} ms -> {ffcx_us_per_call:.3f} us/call",
        flush=True,
    )
    ffcx_total_s = ffcx_compile_s + ffcx_calls_s
    print(f"compile + {n_calls} calls, cold start: {ffcx_total_s * 1e3:.3f} ms", flush=True)
    print(flush=True)

    # rows: list of (row_label, compile_s, us_per_call, total_s)
    rows = [("FFCx (legacy UFL)", ffcx_compile_s, ffcx_us_per_call, ffcx_total_s)]

    gen_label = "UFLx -> MLIR"
    for pv_label, pv_attr in PIPELINE_VARIANTS:
        pipeline = None if pv_attr is None else getattr(mlir_harness, pv_attr)
        compile_s, us_per_call, total_s, best = run_variant(
            f"{gen_label} [{pv_label}]",
            lambda degree, _p=pipeline: uflx_compile(degree, w_dofs, _p),
            degree,
            n_calls,
            vec_ref,
        )
        print(flush=True)
        rows.append((f"{gen_label} [{pv_label}/{best}]", compile_s, us_per_call, total_s))

    print(f"--- Summary (degree={degree}, N={n_calls}) ---", flush=True)
    print(
        f"{'path':<52} {'compile (ms)':>14} {'us/call':>12} {'cold total (ms)':>18}",
        flush=True,
    )
    for label, compile_s, us_per_call, total_s in rows:
        print(
            f"{label:<52} {compile_s * 1e3:>14.3f} {us_per_call:>12.3f} {total_s * 1e3:>18.3f}",
            flush=True,
        )
    print(flush=True)

    for label, compile_s, us_per_call, total_s in rows[1:]:
        print(
            f"FFCx / [{label}] -- compile: {ffcx_compile_s / compile_s:.3f}x  "
            f"per-call: {ffcx_us_per_call / us_per_call:.3f}x  "
            f"cold total: {ffcx_total_s / total_s:.3f}x",
            flush=True,
        )
    print(flush=True)
    print(
        "Linear action of inner(grad(w), grad(v))*dx for a Coefficient w -- the local", flush=True
    )
    print(
        "tensor is a length-ndofs vector (stiffness_matrix @ w's dof vector), not a matrix,",
        flush=True,
    )
    print(
        "and both generated kernels take a third coefficients-buffer argument alongside the",
        flush=True,
    )
    print(
        "local tensor and cell coordinates. Otherwise identical methodology to",
        flush=True,
    )
    print(
        "ffcx_compare_uflx.py: compile time covers the WHOLE UFLx side of the pipeline (form ->",
        flush=True,
    )
    print(
        "quadrature/geometry/tabulation lowering -> MLIR generation -> JIT) against FFCx's",
        flush=True,
    )
    print(
        "compile_forms() (UFL analysis + C codegen + C compile); per-call numbers reuse",
        flush=True,
    )
    print(
        "pre-built buffers/descriptors on both sides throughout. [baseline] uses UFLX_PIPELINE",
        flush=True,
    )
    print(
        "(plain dialect conversion); [licm] additionally runs canonicalize/cse/",
        flush=True,
    )
    print(
        "loop-invariant-code-motion first -- compare the two rows to see whether that",
        flush=True,
    )
    print("closes any remaining gap.", flush=True)


if __name__ == "__main__":
    main()
