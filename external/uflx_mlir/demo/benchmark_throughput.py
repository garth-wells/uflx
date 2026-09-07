"""Benchmark demo/assemble_mesh_gpu.py's single-call assembly kernel's
throughput (dofs/sec) across a log-spaced range of mesh sizes, and plot
the results.

Reuses assemble_mesh_gpu.py's own build_mesh/build_dofmap/
assemble_global_matrix[_gpu] rather than reimplementing any of it --
imported the same way that module imports its own sibling `harness`
(sys.path.insert(0, this directory), then a plain `import`, since demo/
is a loose collection of scripts, not a package).

Usage:
    python3 demo/benchmark_throughput.py [--degree N] [--backend cpu|gpu|both]
        [--cubin-chip CHIP] [--n-min N] [--n-max N] [--num-points N]
        [--csv PATH] [--plot PATH]

    Mesh sizes are log-spaced integers from --n-min to --n-max inclusive
    (default 5 to 60 -- large enough that the largest few points are well
    past 1e6 global dofs, per the point of this benchmark: showing
    throughput once GPU launch latency/CPU startup overhead is no longer
    the dominant cost). Duplicate roundings collapse (e.g. small --n-min
    with many --num-points), so the actual point count can come out below
    --num-points.

    Needs a working uflx/MLIR install -- this can only run wherever that's
    available (this repo's own environment), same as assemble_mesh_gpu.py
    itself.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import assemble_mesh_gpu as amg  # noqa: E402  (see sys.path.insert above)


def log_spaced_mesh_sizes(n_min: int, n_max: int, num_points: int) -> list[int]:
    """Distinct integer mesh resolutions n, log-spaced from n_min to
    n_max inclusive. May return fewer than num_points if rounding to the
    nearest integer collapses two requested points onto the same n (more
    likely near n_min, where geomspace's points are closest together)."""
    if n_min <= 0 or n_max <= 0 or n_min > n_max:
        raise ValueError(f"need 0 < n_min <= n_max, got n_min={n_min}, n_max={n_max}")
    raw = np.geomspace(n_min, n_max, num_points)
    return sorted({int(round(x)) for x in raw})


def run_backend(
    n_values: list[int], degree: int, backend: str, cubin_chip: str
) -> list[dict]:
    """Run one backend across every mesh size in n_values, returning one
    result dict per size. Mirrors main()'s own per-size sequence (build
    mesh, build dofmap, assemble) but keeps mesh/dofmap build timing
    separate from the single kernel call/launch timing, since only the
    latter is what "dofs/sec assembled" (the user's own framing) means --
    see assemble_global_matrix's/assemble_global_matrix_gpu's own
    return_timing option, added alongside this script specifically so
    this doesn't have to reimplement or parse-back-out either kernel
    call's own timing.
    """
    if backend == "gpu":
        assemble_fn = amg.assemble_global_matrix_gpu
        assemble_kwargs = {"cubin_chip": cubin_chip, "return_timing": True}
    else:
        assemble_fn = amg.assemble_global_matrix
        assemble_kwargs = {"return_timing": True}

    results = []
    for n in n_values:
        ncells = 6 * n**3
        t_pre0 = time.perf_counter()
        coords, cells = amg.build_mesh(n)
        cell_dofs, ndofs, ndofs_global = amg.build_dofmap(cells, len(coords), degree)
        t_pre1 = time.perf_counter()

        avals, acols, arowptr, elapsed = assemble_fn(
            coords, cells, cell_dofs, ndofs, ndofs_global, degree, **assemble_kwargs
        )
        dofs_per_sec = amg._dofs_per_sec(ndofs_global, elapsed)

        result = dict(
            backend=backend,
            n=n,
            ncells=ncells,
            ndofs_global=ndofs_global,
            nnz=len(acols),
            preprocessing_seconds=t_pre1 - t_pre0,
            kernel_seconds=elapsed,
            dofs_per_sec=dofs_per_sec,
        )
        results.append(result)
        print(
            f"[{backend}] n={n}: {ncells} cells, {ndofs_global} dofs, "
            f"{elapsed:.4f}s, {dofs_per_sec:.3e} dofs/sec "
            f"(mesh+dofmap build: {t_pre1 - t_pre0:.2f}s)"
        )
    return results


def make_plot(results: list[dict], plot_path: str) -> None:
    """Two side-by-side log-log panels: throughput (the user's own ask)
    and raw kernel wall-clock time (context for it), both vs global dof
    count -- x axis is ndofs_global rather than mesh resolution n, since
    that's the quantity the user's own "more than 1e6 dofs" framing and
    _dofs_per_sec are both defined in terms of. A vertical marker at 1e6
    dofs flags the threshold the user gave for GPU launch latency (or,
    on the CPU backend, Python/MLIR call overhead) no longer being a
    factor."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_throughput, ax_time) = plt.subplots(1, 2, figsize=(13, 5.5))

    backends = sorted({r["backend"] for r in results})
    colors = {"cpu": "tab:blue", "gpu": "tab:red"}
    for backend in backends:
        rows = sorted(
            (r for r in results if r["backend"] == backend),
            key=lambda r: r["ndofs_global"],
        )
        ndofs = [r["ndofs_global"] for r in rows]
        dofs_per_sec = [r["dofs_per_sec"] for r in rows]
        kernel_seconds = [r["kernel_seconds"] for r in rows]
        color = colors.get(backend)
        ax_throughput.plot(
            ndofs, dofs_per_sec, marker="o", label=backend.upper(), color=color
        )
        ax_time.plot(
            ndofs, kernel_seconds, marker="o", label=backend.upper(), color=color
        )

    for ax in (ax_throughput, ax_time):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Global dofs")
        ax.axvline(1e6, color="gray", linestyle=":", linewidth=1)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend()

    ax_throughput.set_ylabel("Throughput (dofs/sec)")
    ax_throughput.set_title("Assembly throughput vs mesh size")
    ax_throughput.text(
        1e6,
        ax_throughput.get_ylim()[0],
        " 10⁶ dofs",
        color="gray",
        fontsize=8,
        va="bottom",
    )

    ax_time.set_ylabel("Single kernel call/launch time (s)")
    ax_time.set_title("Assembly wall-clock time vs mesh size")

    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--backend", choices=["cpu", "gpu", "both"], default="cpu")
    parser.add_argument("--cubin-chip", default="sm_80")
    parser.add_argument("--n-min", type=int, default=5)
    parser.add_argument("--n-max", type=int, default=60)
    parser.add_argument("--num-points", type=int, default=12)
    parser.add_argument("--csv", default="throughput.csv")
    parser.add_argument("--plot", default="throughput.png")
    args = parser.parse_args()

    if args.backend in ("gpu", "both"):
        cuda_runtime_lib = amg._find_cuda_runtime_lib()
        if not cuda_runtime_lib:
            raise SystemExit(
                "GPU backend requested but libmlir_cuda_runtime.so was not "
                "found -- rebuild MLIR with -DMLIR_ENABLE_CUDA_RUNNER=ON, or "
                "set MLIR_CUDA_RUNTIME_LIB to its path."
            )

    n_values = log_spaced_mesh_sizes(args.n_min, args.n_max, args.num_points)
    print(f"mesh sizes (n): {n_values}")
    print(
        f"largest ({n_values[-1]}): "
        f"{6 * n_values[-1] ** 3} cells -- see this script's own docstring "
        f"for why sizes are log-spaced and what the plotted x axis is.\n"
    )

    backends = ["cpu", "gpu"] if args.backend == "both" else [args.backend]
    all_results: list[dict] = []
    for backend in backends:
        all_results.extend(run_backend(n_values, args.degree, backend, args.cubin_chip))

    with open(args.csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_results[0].keys()))
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nwrote {args.csv}")

    make_plot(all_results, args.plot)
    print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
