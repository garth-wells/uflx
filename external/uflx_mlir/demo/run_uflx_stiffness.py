"""UFLx form -> MLIR (uflx_mlir.emit.generate_mlir_module, Python op-builder
API) -> JIT'd P{degree} tetrahedron Laplacian-stiffness kernel, for any
Lagrange degree >= 1 -- not just P1.

This is the first end-to-end proof that "UFLx form -> MLIR kernel" works,
reusing:
  - uflx_codegeneration's existing quadrature/geometry/tabulation pipeline
    (unchanged, just driven with a tetrahedron quadrature rule instead of
    its hardcoded triangle one)
  - harness.py's existing parse -> lower -> JIT -> packed-calling-convention
    machinery (unchanged)
  - generate_kernel.reference_stiffness as the independent correctness
    check, for every degree (that function is itself cross-checked, in
    generate_kernel.main(), against the closed-form P1 reference before
    being trusted for degree>=2).

uflx_mlir is no longer a local prototype living in this repo -- it's the
sibling package ../uflx_mlir (see ../uflx_mlir/emit.py), installed as part
of this same uflx checkout. That package builds MLIR directly via the
Python op-builder API (mlir.ir.Operation.create); this script builds the
ExecutionEngine straight from the in-memory Module via
harness.build_engine_from_module(), with no MLIR-text round-trip at all.
(This script used to have a sibling, run_uflx_p1_ops.py, covering exactly
this op-builder path while this file drove a separate text-emitting
generator and JIT'd the text via a write/re-parse round-trip -- the two
generators were consolidated upstream into the single generate_mlir_module()
this repo now uses, the text round-trip stopped pulling its weight once the
op-builder path was validated end-to-end, and the sibling script was
folded into this one.)

Neither uflx_mlir.emit.generate_mlir_module nor the pipeline it drives
hardcodes a DOF count or loop bound anywhere -- ndofs and every constant
needed as a loop bound/index are read off the lowered graph
(AddToLocalTensor.shape, collect_int_constants()) at generation time. Only
this driver hardcodes degree=1 by default (element degree, %A's 4x4 shape,
and the P1-only reference), so that's all that changes here. Geometry
stays a P1 (affine) coordinate map regardless of solution degree, same as
generate_kernel.py's higher-order kernels.

Run:
    python3 demo/run_uflx_stiffness.py        # degree 1 (default)
    python3 demo/run_uflx_stiffness.py 3      # degree 3

Requires uflx_mlir installed (it's the sibling package one level up --
see the package README's "Building LLVM/MLIR with Python bindings"
section for the one-time setup):
    pip install -e ~/Work/uflx/external/uflx_mlir
"""

import sys
from pathlib import Path

import numpy as np

import basix
from basix_uflx import element
from uflx import coordinate_element, dx, function_space, TestFunction, TrialFunction, grad, inner  # noqa: F401

sys.path.insert(0, str(Path(__file__).parent))
import harness as mlir_harness
import generate_kernel
from uflx_mlir.emit import generate_mlir_module


def build_stiffness_form(degree: int):
    """inner(grad(u), grad(v)) * dx on a P{degree} Lagrange tetrahedron
    space. The coordinate map is always P1 (affine), independent of the
    solution degree -- see module docstring. Returns (form, ndofs)."""
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    u = TrialFunction(space)
    v = TestFunction(space)
    return inner(grad(u), grad(v)) * dx, e.dim


def main():
    try:
        degree = int(sys.argv[1]) if len(sys.argv) == 2 else 1
    except ValueError:
        degree = -1
    if degree < 1:
        print("usage: python3 run_uflx_stiffness.py <degree>  (degree >= 1, default 1)")
        sys.exit(1)

    kernel_name = f"tabulate_tensor_p{degree}_stiffness_uflx"
    form, ndofs = build_stiffness_form(degree)

    print(f"Lowering UFLx P{degree} form ({ndofs} dofs) and building MLIR via the op-builder API...")
    module = generate_mlir_module(
        form, degree=degree, kernel_name=kernel_name, cell=basix.CellType.tetrahedron, inline_geometry=True
    )
    print("OK -- module built and verified (module.operation.verify() passed)")

    out_path = Path(__file__).parent / "kernels" / f"p{degree}_stiffness_uflx.mlir"
    out_path.write_text(str(module))
    print(f"wrote {out_path}")

    print("Building ExecutionEngine directly from the in-memory module...")
    engine = mlir_harness.build_engine_from_module(module, mlir_harness.UFLX_PIPELINE)
    caller = mlir_harness.build_caller_for(engine, kernel_name)

    coords = np.array(
        [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
        dtype=np.float64,
    )
    A = np.zeros((ndofs, ndofs), dtype=np.float64)
    caller(A, coords)

    A_ref = generate_kernel.reference_stiffness(coords, degree)

    print(f"UFLx-generated MLIR result (P{degree}, {ndofs}x{ndofs}):\n", A)
    print("Independent basix-quadrature reference result:\n", A_ref)
    np.testing.assert_allclose(A, A_ref, rtol=1e-9, atol=1e-8)
    print(f"MATCH: UFLx-driven MLIR kernel agrees with the independent P{degree} reference")

    row_sums = A.sum(axis=1)
    print("Row sums (should be ~0):", row_sums)
    np.testing.assert_allclose(row_sums, 0.0, atol=1e-8)
    print("Row-sum check OK")


if __name__ == "__main__":
    main()
