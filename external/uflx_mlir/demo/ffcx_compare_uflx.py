"""Full-pipeline comparison: UFLx form -> MLIR -> JIT'd kernel, versus
legacy-UFL form -> FFCx -> compiled kernel, for a CG-`degree` Laplacian
stiffness kernel on a tetrahedron.

Unlike ffcx_compare.py (which times JIT-ing an ALREADY-GENERATED .mlir
kernel file, from generate_kernel.py's basix-driven quadrature loop, against
FFCx), this measures the WHOLE UFLx
side as "compile" time too: UFLx form -> uflx_codegeneration's
quadrature/geometry/tabulation lowering -> MLIR generation -> JIT.
uflx_mlir.emit.generate_mlir_module (the sibling package ../uflx_mlir,
installed as part of this same uflx checkout) builds the module via the
Python op-builder API and uflx_compile() JITs it directly from the
in-memory Module, with no
MLIR-text round-trip -- that used to be one of two paths compared here (a
text-emitting generator was JIT'd via str(module) -> parse, alongside this
op-builder path); the text generator was consolidated away upstream and
the text round-trip stopped pulling its weight once this path was
validated end-to-end, so only this one remains. This is the fair
apples-to-apples comparison against FFCx's compile_forms(),
which does UFL analysis + C codegen + a real C compile in one call --
ffcx_compare.py's "compile" number, by contrast, only measures parsing an
already-generated MLIR kernel, since generate_kernel.py's own generation
step (basix tabulation + string templating) happens separately beforehand
and isn't timed there.

Steady-state per-call timing and the pre-built-buffers/descriptors
methodology (buffers/descriptors built once, reused across all calls) are
unchanged from ffcx_compare.py -- see that file's docstring for the
reasoning.

Run (same venv, fenics-ffcx installed):
    python3 demo/ffcx_compare_uflx.py        # degree 1 (default)
    python3 demo/ffcx_compare_uflx.py 2
    python3 demo/ffcx_compare_uflx.py 3 50000  # override N_CALLS
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
from run_uflx_stiffness import build_stiffness_form

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


def build_legacy_ufl_form(degree: int):
    """The same form as run_uflx_stiffness.build_stiffness_form (P{degree} Lagrange
    stiffness, P1 affine coordinate map), but expressed in legacy UFL --
    what FFCx actually consumes. Must use the same Lagrange variant
    (equispaced) as the UFLx side -- see ffcx_compare.py's build_form()
    docstring for why this matters starting at degree 3.
    """
    element = basix.ufl.element(
        "Lagrange",
        "tetrahedron",
        degree,
        lagrange_variant=basix.LagrangeVariant.equispaced,
    )
    coord_element = basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,))
    domain = ufl.Mesh(coord_element)
    V = ufl.FunctionSpace(domain, element)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    return ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx


def ffcx_compile(degree: int):
    """Returns (compile_seconds, fast_call, A). Functionally identical to
    ffcx_compare.py's ffcx_compile() -- kept self-contained here rather
    than imported, so this script doesn't couple to that one's differently
    -scoped build_form().
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
    A = np.zeros((ndofs, ndofs), dtype=np.float64)

    A_buf = ffi.from_buffer("double[]", A)
    coords_buf = ffi.from_buffer("double[]", coords_flat)
    eli_buf = ffi.from_buffer("int[]", entity_local_index)
    qp_buf = ffi.from_buffer("uint8_t[]", quadrature_permutation)

    def fast_call():
        tabulate_tensor(A_buf, ffi.NULL, ffi.NULL, coords_buf, eli_buf, qp_buf, ffi.NULL)

    shutil.rmtree(cache_dir, ignore_errors=True)
    return t1 - t0, fast_call, A


def _mlir_fast_calls(engine, kernel_name: str, ndofs: int):
    """Same dual-variant ('invoke' / 'direct_ctypes') call-builder as
    ffcx_compare.py's mlir_compile() -- factored out here since both UFLx
    variants below need it identically.
    """
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
        print(
            f"(direct ctypes lookup path unavailable: {e!r} -- skipping that variant)", flush=True
        )

    return fast_calls, A


def uflx_compile(degree: int, pipeline: str | None = None):
    """Full pipeline, direct-from-module: UFLx form -> generate_mlir_module's
    in-memory Module -> JIT, with no MLIR-text round-trip at all.
    compile_seconds covers all of that -- the UFLx-side equivalent of
    FFCx's compile_forms() timing.

    pipeline defaults to mlir_harness.UFLX_PIPELINE (plain dialect
    conversion, no optimization passes); pass
    mlir_harness.UFLX_OPTIMIZED_PIPELINE to additionally run
    canonicalize/cse/loop-invariant-code-motion first -- see that
    constant's docstring for why.
    """
    pipeline = pipeline if pipeline is not None else mlir_harness.UFLX_PIPELINE
    form, ndofs = build_stiffness_form(degree)
    kernel_name = f"tabulate_tensor_p{degree}_stiffness_uflx"

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

    fast_calls, A = _mlir_fast_calls(engine, kernel_name, ndofs)
    return t1 - t0, fast_calls, A


def time_calls(fast_call, n):
    t0 = time.perf_counter()
    for _ in range(n):
        fast_call()
    t1 = time.perf_counter()
    return t1 - t0


def run_variant(label: str, compile_fn, degree: int, n_calls: int, A_ref: np.ndarray):
    """Runs one 'compile, validate, time steady-state calls' pass -- shared
    by FFCx and the UFLx compile path below. Returns
    (compile_seconds, best_us_per_call, best_total_seconds, best_variant_name).
    """
    print(f"--- {label} ---", flush=True)
    compile_s, fast_calls, A = compile_fn(degree)
    print(f"compile: {compile_s * 1e3:.3f} ms", flush=True)

    results = {}
    for variant, fast_call in fast_calls.items():
        # A is shared across every variant in fast_calls (see
        # _mlir_fast_calls) -- now that the UFLx kernel no longer
        # zero-initializes itself (it's pure accumulate, matching FFCx's/
        # UFC's tabulate_tensor convention), a variant that runs after
        # another has already accumulated into A would otherwise fail this
        # check against a stale, non-zero buffer. Reset here, once per
        # variant, entirely outside the timed region below -- time_calls
        # itself still reuses A across all n_calls without rezeroing,
        # unchanged from before.
        A.fill(0.0)
        fast_call()
        np.testing.assert_allclose(A, A_ref, rtol=1e-9, atol=1e-8)
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


# (label, pipeline) pairs run for the UFLx compile path below. "baseline" is
# plain dialect conversion, no optimization; "licm" additionally runs
# canonicalize/cse/loop-invariant-code-motion first -- see
# mlir_harness.UFLX_OPTIMIZED_PIPELINE's docstring for why this was added
# (a ~27x per-call slowdown vs FFCx at degree 3, traced to the generated
# loop nest -- i/j/quadrature-point, quadrature innermost -- flattening the
# WHOLE per-entry expression, including the cell's Jacobian/detJ/cofactor
# terms, into the innermost loop body with no hoisting of the parts that
# don't depend on the inner loop variables).
PIPELINE_VARIANTS = [
    ("baseline", None),  # None -> uflx_compile's own UFLX_PIPELINE default
    ("licm", "UFLX_OPTIMIZED_PIPELINE"),  # resolved via getattr below
]


def main():
    degree = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n_calls = int(sys.argv[2]) if len(sys.argv) > 2 else N_CALLS_DEFAULT
    print(f"=== degree {degree}, N_CALLS={n_calls} ===", flush=True)
    A_ref = generate_kernel.reference_stiffness(COORDS, degree)

    print("--- FFCx (legacy UFL) ---", flush=True)
    ffcx_compile_s, ffcx_call, A_ffcx = ffcx_compile(degree)
    ffcx_call()
    np.testing.assert_allclose(A_ffcx, A_ref, rtol=1e-9, atol=1e-8)
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
            lambda degree, _p=pipeline: uflx_compile(degree, _p),
            degree,
            n_calls,
            A_ref,
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
    print("compile time here covers the WHOLE UFLx side of the pipeline (form ->", flush=True)
    print(
        "quadrature/geometry/tabulation lowering -> MLIR generation -> JIT) -- the fair", flush=True
    )
    print(
        "comparison against FFCx's compile_forms() (UFL analysis + C codegen + C compile).",
        flush=True,
    )
    print(
        "Per-call numbers reuse pre-built buffers/descriptors on both sides throughout,", flush=True
    )
    print(
        "same methodology as ffcx_compare.py. [baseline] uses UFLX_PIPELINE (plain dialect",
        flush=True,
    )
    print(
        "conversion); [licm] additionally runs canonicalize/cse/loop-invariant-code-motion",
        flush=True,
    )
    print("first -- compare the two rows to see whether that closes the gap.", flush=True)


if __name__ == "__main__":
    main()
