"""Generate a P{degree} Lagrange Laplacian-stiffness kernel (a real
quadrature-point loop, using scf.for) as MLIR text, from basix's own
reference-element tabulation -- the same quadrature points/weights and
basis-function-gradient tables FFCx's own code generator would use for this
form, on a tetrahedron.

Run (same venv as harness.py, after `pip install fenics-ffcx` for basix):
    python3 demo/generate_kernel.py 1   # writes kernels/p1_stiffness.mlir
    python3 demo/generate_kernel.py 2   # writes kernels/p2_stiffness.mlir
    python3 demo/generate_kernel.py 3   # writes kernels/p3_stiffness.mlir

Geometry is a P1-mapped (affine) tetrahedron regardless of the solution
element's degree. The generated P1 kernel uses the same general quadrature-loop
structure as higher degrees; its basis gradients are constant and the loop has
one quadrature point. The simplifying assumption is positively-oriented
tetrahedra (uses detJ directly, not abs(detJ)).

Uses the equispaced Lagrange variant (textbook nodal points), not DOLFINx's
default (gll_isaac) -- simpler to reason about; swap LAGRANGE_VARIANT below
for exact parity with a real dolfinx.fem.functionspace. ffcx_compare.py's
build_form() must use the SAME variant -- for degree<=2 any variant coincides
(<=1 point per edge), but degree>=3 has multiple points per edge, where
equispaced and other variants genuinely differ (found the hard way: this is
exactly what broke the P3 triangle comparison before the 3D rework).
"""

import sys
from pathlib import Path

import basix
import numpy as np

LAGRANGE_VARIANT = basix.LagrangeVariant.equispaced
CELL = basix.CellType.tetrahedron

KERNELS_DIR = Path(__file__).parent / "kernels"


def _make_quadrature(cell, degree):
    """basix.make_quadrature's signature has changed across versions -- try
    the plain (cell, degree) form first, fall back to the (type, cell, degree)
    form some versions require.
    """
    try:
        return basix.make_quadrature(cell, degree)
    except TypeError:
        return basix.make_quadrature(basix.QuadratureType.default, cell, degree)


def tabulate(degree: int):
    """Returns (weights[nq], dphi_dx[nq, ndofs], dphi_dy[nq, ndofs],
    dphi_dz[nq, ndofs]) for the degree-`degree` Lagrange element on a
    reference tetrahedron, with defensive checks on basix's tabulate()
    output shape/derivative-ordering, since that's exactly the kind of thing
    that's easy to get subtly wrong and this script has no way to be
    test-run before you run it.
    """
    element = basix.create_element(basix.ElementFamily.P, CELL, degree, LAGRANGE_VARIANT)
    ndofs = element.dim

    # grad(phi_i).grad(phi_j) has degree 2*(degree-1) -- quadrature must be
    # exact to at least that.
    qdeg = max(2 * (degree - 1), 1)
    points, weights = _make_quadrature(CELL, qdeg)
    nq = len(weights)

    # tab[deriv_index, point, dof, component]; for a 3D cell with nderivs=1,
    # deriv_index 0/1/2/3 = value/d_dx/d_dy/d_dz -- this ordering is asserted
    # below, not just assumed.
    tab = element.tabulate(1, points)
    values = tab[0, :, :, 0]
    dphi_dx = tab[1, :, :, 0]
    dphi_dy = tab[2, :, :, 0]
    dphi_dz = tab[3, :, :, 0]

    assert values.shape == (nq, ndofs), f"unexpected tabulate() shape {values.shape}"

    # Partition-of-unity invariants: sum_i phi_i == 1 everywhere, so its
    # derivatives are 0 everywhere. Cheap, strong checks on the deriv-index
    # ordering assumption above.
    assert np.allclose(values.sum(axis=1), 1.0), (
        "basis values don't sum to 1 at quadrature points -- "
        "tabulate() indexing assumption is wrong"
    )
    for name, d in (("dx", dphi_dx), ("dy", dphi_dy), ("dz", dphi_dz)):
        assert np.allclose(d.sum(axis=1), 0.0, atol=1e-10), (
            f"d/{name} of the basis sum isn't 0 -- tabulate() deriv-index "
            f"ordering assumption is wrong"
        )

    return weights, dphi_dx, dphi_dy, dphi_dz


def reference_stiffness(coords: np.ndarray, degree: int) -> np.ndarray:
    """General quadrature-based numpy/basix reference for CG-`degree`
    Laplacian stiffness on an affine-mapped tetrahedron -- independent of
    the MLIR kernel (uses numpy's own 3x3 inverse/det, not the hand-derived
    cofactor formulas the generated kernel uses), so it can validate it.
    For degree=1 this must reduce to exactly the same answer as
    harness.reference_p1_stiffness; that cross-check runs in __main__ below
    before this function is trusted for degree>=2.
    """
    weights, dphi_dx, dphi_dy, dphi_dz = tabulate(degree)

    x0, x1, x2, x3 = coords
    J = np.column_stack([x1 - x0, x2 - x0, x3 - x0])  # 3x3
    detJ = np.linalg.det(J)
    Jinv = np.linalg.inv(J)

    ndofs = dphi_dx.shape[1]
    A = np.zeros((ndofs, ndofs))
    for q, w in enumerate(weights):
        ref_grads = np.stack([dphi_dx[q], dphi_dy[q], dphi_dz[q]], axis=1)  # [dof, 3]
        phys_grads = ref_grads @ Jinv
        A += w * detJ * (phys_grads @ phys_grads.T)
    return A


def _fmt1(arr) -> str:
    return "[" + ", ".join(repr(float(v)) for v in arr) + "]"


def _fmt2(arr) -> str:
    return "[" + ", ".join(_fmt1(row) for row in arr) + "]"


TEMPLATE = """\
// P@@DEGREE@@ (degree @@DEGREE@@) Laplacian stiffness matrix on a single
// tetrahedron, affine map -- generated by demo/generate_kernel.py from
// basix's own reference-element tabulation (quadrature points/weights,
// basis-function gradients), the same numerical inputs FFCx's own code
// generator uses for this form.
//
// The same quadrature-loop structure is used at every degree. For P1, the
// basis gradients are constant and the loop contains one quadrature point;
// at higher degrees the gradients vary over the cell.
//
// Assumes positively-oriented tetrahedra (uses detJ directly, not
// abs(detJ)) -- same simplifying assumption as the P1 kernel. Lagrange
// variant: equispaced.

memref.global "private" constant @p@@DEGREE@@_weights : memref<@@NQ@@xf64> = dense<@@WEIGHTS@@>
memref.global "private" constant @p@@DEGREE@@_dphidx_ref : memref<@@NQ@@x@@NDOFS@@xf64> = dense<@@DPHIDX@@>
memref.global "private" constant @p@@DEGREE@@_dphidy_ref : memref<@@NQ@@x@@NDOFS@@xf64> = dense<@@DPHIDY@@>
memref.global "private" constant @p@@DEGREE@@_dphidz_ref : memref<@@NQ@@x@@NDOFS@@xf64> = dense<@@DPHIDZ@@>

func.func @@@KNAME@@(%A: memref<@@NDOFS@@x@@NDOFS@@xf64>, %coords: memref<4x3xf64>)
    attributes { llvm.emit_c_interface } {
  %c0 = arith.constant 0 : index
  %c1 = arith.constant 1 : index
  %c2 = arith.constant 2 : index
  %c3 = arith.constant 3 : index
  %cN = arith.constant @@NDOFS@@ : index
  %cNQ = arith.constant @@NQ@@ : index
  %zero = arith.constant 0.0 : f64

  %x0 = memref.load %coords[%c0, %c0] : memref<4x3xf64>
  %y0 = memref.load %coords[%c0, %c1] : memref<4x3xf64>
  %z0 = memref.load %coords[%c0, %c2] : memref<4x3xf64>
  %x1 = memref.load %coords[%c1, %c0] : memref<4x3xf64>
  %y1 = memref.load %coords[%c1, %c1] : memref<4x3xf64>
  %z1 = memref.load %coords[%c1, %c2] : memref<4x3xf64>
  %x2 = memref.load %coords[%c2, %c0] : memref<4x3xf64>
  %y2 = memref.load %coords[%c2, %c1] : memref<4x3xf64>
  %z2 = memref.load %coords[%c2, %c2] : memref<4x3xf64>
  %x3 = memref.load %coords[%c3, %c0] : memref<4x3xf64>
  %y3 = memref.load %coords[%c3, %c1] : memref<4x3xf64>
  %z3 = memref.load %coords[%c3, %c2] : memref<4x3xf64>

  %j00 = arith.subf %x1, %x0 : f64
  %j01 = arith.subf %x2, %x0 : f64
  %j02 = arith.subf %x3, %x0 : f64
  %j10 = arith.subf %y1, %y0 : f64
  %j11 = arith.subf %y2, %y0 : f64
  %j12 = arith.subf %y3, %y0 : f64
  %j20 = arith.subf %z1, %z0 : f64
  %j21 = arith.subf %z2, %z0 : f64
  %j22 = arith.subf %z3, %z0 : f64

  // cf[a][b] == detJ * Jinv[a][b] (the (a,b) cofactor of J, transposed
  // appropriately -- see kernels/p1_stiffness.mlir's derivation comment for
  // the full standard 3x3-inverse formula this follows).
  %cf00_a = arith.mulf %j11, %j22 : f64
  %cf00_b = arith.mulf %j12, %j21 : f64
  %cf00 = arith.subf %cf00_a, %cf00_b : f64

  %cf01_a = arith.mulf %j02, %j21 : f64
  %cf01_b = arith.mulf %j01, %j22 : f64
  %cf01 = arith.subf %cf01_a, %cf01_b : f64

  %cf02_a = arith.mulf %j01, %j12 : f64
  %cf02_b = arith.mulf %j02, %j11 : f64
  %cf02 = arith.subf %cf02_a, %cf02_b : f64

  %cf10_a = arith.mulf %j12, %j20 : f64
  %cf10_b = arith.mulf %j10, %j22 : f64
  %cf10 = arith.subf %cf10_a, %cf10_b : f64

  %cf11_a = arith.mulf %j00, %j22 : f64
  %cf11_b = arith.mulf %j02, %j20 : f64
  %cf11 = arith.subf %cf11_a, %cf11_b : f64

  %cf12_a = arith.mulf %j02, %j10 : f64
  %cf12_b = arith.mulf %j00, %j12 : f64
  %cf12 = arith.subf %cf12_a, %cf12_b : f64

  %cf20_a = arith.mulf %j10, %j21 : f64
  %cf20_b = arith.mulf %j11, %j20 : f64
  %cf20 = arith.subf %cf20_a, %cf20_b : f64

  %cf21_a = arith.mulf %j01, %j20 : f64
  %cf21_b = arith.mulf %j00, %j21 : f64
  %cf21 = arith.subf %cf21_a, %cf21_b : f64

  %cf22_a = arith.mulf %j00, %j11 : f64
  %cf22_b = arith.mulf %j01, %j10 : f64
  %cf22 = arith.subf %cf22_a, %cf22_b : f64

  // detJ = j00*cf00 + j01*cf10 + j02*cf20 (cofactor expansion along row 0)
  %detJ_t0 = arith.mulf %j00, %cf00 : f64
  %detJ_t1 = arith.mulf %j01, %cf10 : f64
  %detJ_t2 = arith.mulf %j02, %cf20 : f64
  %detJ_s0 = arith.addf %detJ_t0, %detJ_t1 : f64
  %detJ = arith.addf %detJ_s0, %detJ_t2 : f64

  %weights = memref.get_global @p@@DEGREE@@_weights : memref<@@NQ@@xf64>
  %dphidx = memref.get_global @p@@DEGREE@@_dphidx_ref : memref<@@NQ@@x@@NDOFS@@xf64>
  %dphidy = memref.get_global @p@@DEGREE@@_dphidy_ref : memref<@@NQ@@x@@NDOFS@@xf64>
  %dphidz = memref.get_global @p@@DEGREE@@_dphidz_ref : memref<@@NQ@@x@@NDOFS@@xf64>

  // Scratch space for the per-quadrature-point "numerator" vectors
  // (detJ * physical gradient), indexed by dof. Filled once per q and read
  // by BOTH the i-loop and j-loop -- avoids recomputing num[j] from scratch
  // inside the j-loop for every i (an O(ndofs) redundant-work factor that
  // only mattered once ndofs got big enough for the arithmetic to stop
  // being negligible next to call overhead).
  %numx_scratch = memref.alloca() : memref<@@NDOFS@@xf64>
  %numy_scratch = memref.alloca() : memref<@@NDOFS@@xf64>
  %numz_scratch = memref.alloca() : memref<@@NDOFS@@xf64>

  scf.for %i0 = %c0 to %cN step %c1 {
    scf.for %j0 = %c0 to %cN step %c1 {
      memref.store %zero, %A[%i0, %j0] : memref<@@NDOFS@@x@@NDOFS@@xf64>
    }
  }

  scf.for %q = %c0 to %cNQ step %c1 {
    %w = memref.load %weights[%q] : memref<@@NQ@@xf64>

    scf.for %k = %c0 to %cN step %c1 {
      %dxk = memref.load %dphidx[%q, %k] : memref<@@NQ@@x@@NDOFS@@xf64>
      %dyk = memref.load %dphidy[%q, %k] : memref<@@NQ@@x@@NDOFS@@xf64>
      %dzk = memref.load %dphidz[%q, %k] : memref<@@NQ@@x@@NDOFS@@xf64>

      %numxk_0 = arith.mulf %dxk, %cf00 : f64
      %numxk_1 = arith.mulf %dyk, %cf10 : f64
      %numxk_2 = arith.mulf %dzk, %cf20 : f64
      %numxk_s = arith.addf %numxk_0, %numxk_1 : f64
      %numxk = arith.addf %numxk_s, %numxk_2 : f64

      %numyk_0 = arith.mulf %dxk, %cf01 : f64
      %numyk_1 = arith.mulf %dyk, %cf11 : f64
      %numyk_2 = arith.mulf %dzk, %cf21 : f64
      %numyk_s = arith.addf %numyk_0, %numyk_1 : f64
      %numyk = arith.addf %numyk_s, %numyk_2 : f64

      %numzk_0 = arith.mulf %dxk, %cf02 : f64
      %numzk_1 = arith.mulf %dyk, %cf12 : f64
      %numzk_2 = arith.mulf %dzk, %cf22 : f64
      %numzk_s = arith.addf %numzk_0, %numzk_1 : f64
      %numzk = arith.addf %numzk_s, %numzk_2 : f64

      memref.store %numxk, %numx_scratch[%k] : memref<@@NDOFS@@xf64>
      memref.store %numyk, %numy_scratch[%k] : memref<@@NDOFS@@xf64>
      memref.store %numzk, %numz_scratch[%k] : memref<@@NDOFS@@xf64>
    }

    scf.for %i = %c0 to %cN step %c1 {
      %numxi = memref.load %numx_scratch[%i] : memref<@@NDOFS@@xf64>
      %numyi = memref.load %numy_scratch[%i] : memref<@@NDOFS@@xf64>
      %numzi = memref.load %numz_scratch[%i] : memref<@@NDOFS@@xf64>
      scf.for %j = %c0 to %cN step %c1 {
        %numxj = memref.load %numx_scratch[%j] : memref<@@NDOFS@@xf64>
        %numyj = memref.load %numy_scratch[%j] : memref<@@NDOFS@@xf64>
        %numzj = memref.load %numz_scratch[%j] : memref<@@NDOFS@@xf64>

        %dotx = arith.mulf %numxi, %numxj : f64
        %doty = arith.mulf %numyi, %numyj : f64
        %dotz = arith.mulf %numzi, %numzj : f64
        %dot_s = arith.addf %dotx, %doty : f64
        %dot = arith.addf %dot_s, %dotz : f64
        %wdot = arith.mulf %w, %dot : f64
        %contrib = arith.divf %wdot, %detJ : f64

        %old = memref.load %A[%i, %j] : memref<@@NDOFS@@x@@NDOFS@@xf64>
        %new = arith.addf %old, %contrib : f64
        memref.store %new, %A[%i, %j] : memref<@@NDOFS@@x@@NDOFS@@xf64>
      }
    }
  }
  return
}
"""


def generate_mlir(degree: int) -> str:
    weights, dphi_dx, dphi_dy, dphi_dz = tabulate(degree)
    nq, ndofs = dphi_dx.shape
    text = TEMPLATE
    for token, value in [
        ("@@DEGREE@@", str(degree)),
        ("@@NQ@@", str(nq)),
        ("@@NDOFS@@", str(ndofs)),
        ("@@WEIGHTS@@", _fmt1(weights)),
        ("@@DPHIDX@@", _fmt2(dphi_dx)),
        ("@@DPHIDY@@", _fmt2(dphi_dy)),
        ("@@DPHIDZ@@", _fmt2(dphi_dz)),
        ("@@KNAME@@", f"tabulate_tensor_p{degree}_stiffness"),
    ]:
        text = text.replace(token, value)
    assert "@@" not in text, "leftover unreplaced template token"
    return text


def main():
    try:
        degree = int(sys.argv[1]) if len(sys.argv) == 2 else -1
    except ValueError:
        degree = -1
    if degree < 1:
        print("usage: python3 generate_kernel.py <degree>  (degree >= 1)")
        sys.exit(1)
    if degree > 6:
        print(
            f"note: degree {degree} equispaced Lagrange nodes get increasingly "
            f"ill-conditioned (Runge's phenomenon) -- proceeding anyway, this is a "
            f"numerical-conditioning caveat, not a correctness blocker."
        )

    # Cross-check the general basix-based reference_stiffness() against the
    # already-proven-correct closed-form P1 reference before using it to
    # generate a kernel -- catches a wrong basix API/indexing assumption immediately
    # rather than silently generating a wrong kernel.
    sys.path.insert(0, str(Path(__file__).parent))
    import harness as mlir_harness

    coords = np.array(
        [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
        dtype=np.float64,
    )
    A1_general = reference_stiffness(coords, 1)
    A1_closed_form = mlir_harness.reference_p1_stiffness(coords)
    np.testing.assert_allclose(A1_general, A1_closed_form, rtol=1e-10)
    print("Cross-check OK: general quadrature reference matches the closed-form P1 reference.")

    mlir_text = generate_mlir(degree)
    out_path = KERNELS_DIR / f"p{degree}_stiffness.mlir"
    out_path.write_text(mlir_text)
    print(f"wrote {out_path} ({len(mlir_text.splitlines())} lines)")


if __name__ == "__main__":
    main()
