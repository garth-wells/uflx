"""Benchmark FFCx (legacy UFL) vs UFLx -> MLIR [licm] across CG degrees
1-6 on CPU, and plot compile time and steady-state per-call runtime for
both.

Reuses demo/ffcx_compare_uflx.py's own machinery rather than
reimplementing any of it -- imported the same way that module imports
its own sibling `harness` (sys.path.insert(0, this directory), then a
plain `import`, since demo/ is a loose collection of scripts, not a
package): ffcx_compile() and uflx_compile() do the actual compiling and
build the fast-call closures, run_variant() picks the faster of the
'invoke'/'direct_ctypes' MLIR call variants (see that function's own
docstring), and time_calls() does the steady-state per-call timing --
this script only adds the loop over degree and the plotting.

Only the "licm" pipeline variant is compared here (not "baseline" --
plain dialect conversion with no optimization passes) since
generate_mlir_module's own loop-hoisting (see harness.py's
UFLX_OPTIMIZED_PIPELINE docstring) already closes most of the gap
baseline showed at degree 3+; "licm" is the fairer, faster of the two
UFLx variants to hold up against FFCx.

Usage:
    python3 demo/benchmark_ffcx_vs_mlir.py [--degree-min N] [--degree-max N]
        [--n-calls N] [--csv PATH] [--plot PATH]

    Needs a working uflx/MLIR install AND fenics-ffcx (same venv
    ffcx_compare_uflx.py itself needs) -- this can only run wherever
    that's available (this repo's own environment).
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import ffcx_compare_uflx as fcu  # noqa: E402  (see sys.path.insert above)


def run_degree(degree: int, n_calls: int) -> dict:
    """One degree's worth of ffcx_compare_uflx.py's own main(), minus the
    printing: FFCx compile+validate+time, then UFLx->MLIR [licm]
    compile+validate+time via run_variant (which also picks whichever of
    the 'invoke'/'direct_ctypes' call variants is faster -- see that
    function's own docstring)."""
    print(f"=== degree {degree} ===", flush=True)
    a_ref = fcu.generate_kernel.reference_stiffness(fcu.COORDS, degree)

    ffcx_compile_s, ffcx_call, a_ffcx = fcu.ffcx_compile(degree)
    ffcx_call()
    np.testing.assert_allclose(a_ffcx, a_ref, rtol=1e-9, atol=1e-8)
    ffcx_calls_s = fcu.time_calls(ffcx_call, n_calls)
    ffcx_us_per_call = ffcx_calls_s / n_calls * 1e6
    print(
        f"  FFCx: compile {ffcx_compile_s * 1e3:.3f} ms, "
        f"{ffcx_us_per_call:.3f} us/call",
        flush=True,
    )

    licm_compile_s, licm_us_per_call, _licm_total_s, licm_variant = fcu.run_variant(
        "UFLx -> MLIR [licm]",
        lambda d, _p=fcu.mlir_harness.UFLX_OPTIMIZED_PIPELINE: fcu.uflx_compile(d, _p),
        degree,
        n_calls,
        a_ref,
    )
    print(
        f"  MLIR/licm [{licm_variant}]: compile {licm_compile_s * 1e3:.3f} ms, "
        f"{licm_us_per_call:.3f} us/call",
        flush=True,
    )
    print(flush=True)

    return dict(
        degree=degree,
        ffcx_compile_ms=ffcx_compile_s * 1e3,
        ffcx_us_per_call=ffcx_us_per_call,
        mlir_licm_compile_ms=licm_compile_s * 1e3,
        mlir_licm_us_per_call=licm_us_per_call,
        mlir_licm_variant=licm_variant,
    )


def make_plot(results: list[dict], plot_path: str) -> None:
    """Two side-by-side panels sharing an x axis (CG degree, 1-6):
    compile time (ms) and steady-state per-call runtime (us), each with
    one line for FFCx and one for UFLx -> MLIR [licm]. Y axes are log
    scale -- compile time and runtime both plausibly span more than one
    order of magnitude across degree 1-6 (see ffcx_compare_uflx.py's own
    docstring: a P3 kernel alone showed a ~27x per-call gap before the
    licm/hoisting work)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_compile, ax_runtime) = plt.subplots(1, 2, figsize=(12, 5.5))

    rows = sorted(results, key=lambda r: r["degree"])
    degrees = [r["degree"] for r in rows]

    # Fixed categorical order, not an arbitrary color cycle -- slots 1
    # (blue) and 2 (orange) of the validated default palette, chosen for
    # colorblind-safe adjacent separation (confirmed via this repo's own
    # dataviz skill validator: worst adjacent CVD Delta E 24.7, well
    # clear of the >=8 target) rather than matplotlib's default tab10.
    series = [
        ("FFCx", "ffcx_compile_ms", "ffcx_us_per_call", "#2a78d6"),
        ("MLIR (licm)", "mlir_licm_compile_ms", "mlir_licm_us_per_call", "#eb6834"),
    ]
    for label, compile_key, runtime_key, color in series:
        ax_compile.plot(
            degrees,
            [r[compile_key] for r in rows],
            marker="o",
            markersize=7,
            linewidth=2,
            label=label,
            color=color,
        )
        ax_runtime.plot(
            degrees,
            [r[runtime_key] for r in rows],
            marker="o",
            markersize=7,
            linewidth=2,
            label=label,
            color=color,
        )

    for ax in (ax_compile, ax_runtime):
        ax.set_yscale("log")
        ax.set_xlabel("Lagrange degree")
        ax.set_xticks(degrees)
        ax.grid(True, which="both", ls=":", alpha=0.3, color="#9a9a94")
        ax.legend(frameon=False)

    ax_compile.set_ylabel("Compile time (ms)")
    ax_compile.set_title("Compile time: FFCx vs MLIR (licm)")

    ax_runtime.set_ylabel("Steady-state runtime (µs/call)")
    ax_runtime.set_title("Per-call runtime: FFCx vs MLIR (licm)")

    fig.suptitle("FFCx vs UFLx→MLIR (licm), CPU, degree 1–6")
    fig.tight_layout()
    fig.savefig(plot_path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree-min", type=int, default=1)
    parser.add_argument("--degree-max", type=int, default=6)
    parser.add_argument("--n-calls", type=int, default=fcu.N_CALLS_DEFAULT)
    parser.add_argument("--csv", default="ffcx_vs_mlir.csv")
    parser.add_argument("--plot", default="ffcx_vs_mlir_compare.pdf")
    args = parser.parse_args()

    degrees = list(range(args.degree_min, args.degree_max + 1))
    print(f"degrees: {degrees}, N_CALLS={args.n_calls}\n", flush=True)

    results = [run_degree(degree, args.n_calls) for degree in degrees]

    with open(args.csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    print(f"wrote {args.csv}")

    make_plot(results, args.plot)
    print(f"wrote {args.plot}")

    print("\n--- Summary ---")
    print(
        f"{'degree':>6} {'FFCx compile (ms)':>18} {'MLIR/licm compile (ms)':>24} "
        f"{'FFCx us/call':>14} {'MLIR/licm us/call':>18}"
    )
    for r in results:
        print(
            f"{r['degree']:>6} {r['ffcx_compile_ms']:>18.3f} "
            f"{r['mlir_licm_compile_ms']:>24.3f} {r['ffcx_us_per_call']:>14.3f} "
            f"{r['mlir_licm_us_per_call']:>18.3f}"
        )


if __name__ == "__main__":
    main()
