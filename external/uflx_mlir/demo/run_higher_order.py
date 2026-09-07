"""JIT + validate the quadrature-loop stiffness kernels from generate_kernel.py,
including the generated one-quadrature-point P1 kernel.

Run (after `python3 demo/generate_kernel.py <degree>`):
    python3 demo/run_higher_order.py 1
    python3 demo/run_higher_order.py 2
    python3 demo/run_higher_order.py 3
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import generate_kernel
import harness as mlir_harness


def main():
    try:
        degree = int(sys.argv[1]) if len(sys.argv) == 2 else -1
    except ValueError:
        degree = -1
    if degree < 1:
        print("usage: python3 run_higher_order.py <degree>  (degree >= 1)")
        sys.exit(1)

    kernel_path = Path(__file__).parent / "kernels" / f"p{degree}_stiffness.mlir"
    kernel_name = f"tabulate_tensor_p{degree}_stiffness"
    if not kernel_path.exists():
        print(f"{kernel_path} doesn't exist -- run generate_kernel.py {degree} first")
        sys.exit(1)

    engine = mlir_harness.build_engine_from(kernel_path, mlir_harness.QUADRATURE_PIPELINE)
    caller = mlir_harness.build_caller_for(engine, kernel_name)

    # A scalene tetrahedron, not a reference-aligned one -- less likely to
    # hide an index-swap or transpose bug behind accidental symmetry. Same
    # coordinates as generate_kernel.py's own P1 cross-check.
    coords = np.array(
        [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
        dtype=np.float64,
    )

    _, dphi_dx, _, _ = generate_kernel.tabulate(degree)
    ndofs = dphi_dx.shape[1]
    A = np.zeros((ndofs, ndofs), dtype=np.float64)
    caller(A, coords)

    A_ref = generate_kernel.reference_stiffness(coords, degree)

    print(f"MLIR JIT result (P{degree}, {ndofs}x{ndofs}):\n", A)
    print("Reference result:\n", A_ref)
    np.testing.assert_allclose(
        A, A_ref, rtol=1e-9, atol=1e-8
    )  # atol: near-zero entries shouldn't fail on relative tolerance alone
    print("MATCH: within 1e-9")

    # Independent sanity check (doesn't rely on basix at all): constants are
    # in the kernel of the Laplacian, so every row of a stiffness matrix
    # should sum to ~0. Catches a bug the reference and kernel could share
    # (e.g. both built from the same wrong basix tabulation assumption).
    row_sums = A.sum(axis=1)
    print("Row sums (should be ~0):", row_sums)
    np.testing.assert_allclose(row_sums, 0.0, atol=1e-8)
    print("Row-sum check OK")


if __name__ == "__main__":
    main()
