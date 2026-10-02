"""Check that a generated kernel's MLIR text parses and lowers cleanly through
its pass pipeline, WITHOUT building an ExecutionEngine or JIT'ing anything --
a pure parse+lower step (mlir.ir.Module.parse + mlir.passmanager.PassManager,
what `mlir-opt <file> --pass-pipeline=...` does from the shell, minus the
separate binary). Bad IR or an unknown pass name just raises a Python
exception here; this step was never the crash risk in this repo's history --
that was ExecutionEngine.lookup()'s calling convention, a later, separate step.

Run:
    python3 demo/check_lowering.py 1   # kernels/p1_stiffness.mlir, QUADRATURE_PIPELINE
    python3 demo/check_lowering.py 2   # kernels/p2_stiffness.mlir, QUADRATURE_PIPELINE
    python3 demo/check_lowering.py 3
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import harness as mlir_harness


def main():
    try:
        degree = int(sys.argv[1]) if len(sys.argv) == 2 else -1
    except ValueError:
        degree = -1
    if degree < 1:
        print("usage: python3 check_lowering.py <degree>  (degree >= 1)")
        sys.exit(1)

    kernel_path = Path(__file__).parent / "kernels" / f"p{degree}_stiffness.mlir"
    pipeline = mlir_harness.QUADRATURE_PIPELINE

    print(f"parsing + lowering {kernel_path} with QUADRATURE_PIPELINE ...")
    lowered = mlir_harness.check_lowering(kernel_path, pipeline)
    print(f"OK -- lowered module is {len(lowered.splitlines())} lines. Last 10:")
    print("\n".join(lowered.splitlines()[-10:]))


if __name__ == "__main__":
    main()
