"""Assemble inner(grad(w), grad(v))*dx on CUDA/HIP with automatic cell grouping.

Run from any directory with UFLx and MLIR installed:
    python demo/assemble_linear_gpu.py --degree 1 --n 20 --backend cuda --chip sm_89
    python demo/assemble_linear_gpu.py --degree 2 --n 20 --backend amd --chip gfx1100

The demo's mesh/dof-map builder supports P1/P2. The generator is also tested
at P3/P4. Coefficients are packed per cell from one shared global vector.
Pass --cells-per-block 1 to compare against one cell per block.
"""

from __future__ import annotations

import argparse

import numpy as np
from assemble_mesh_gpu import CELL, _reference_stiffness_cell, build_dofmap, build_mesh
from basix_uflx import element
from uflx import Coefficient, TestFunction, coordinate_element, dx, function_space, grad, inner

from uflx_mlir.gpu_linear import generate_linear_assembly_gpu_module
from uflx_mlir.gpu_runtime import assemble_linear_gpu


def assemble(n: int, degree: int, backend: str, chip: str, cells_per_block: int | None):
    """Assemble a mesh and return its vector and launch measurements."""
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    w, v = Coefficient(space), TestFunction(space)
    form = inner(grad(w), grad(v)) * dx
    module, layout = generate_linear_assembly_gpu_module(
        form, degree, "assemble_linear", CELL, cells_per_block=cells_per_block
    )
    points, cells = build_mesh(n)
    dofs, ndofs, nglobal = build_dofmap(cells, len(points), degree)
    dofs = dofs.reshape(len(cells), ndofs)
    coords = np.ascontiguousarray(points[np.asarray(cells)])
    global_w = np.random.default_rng(812).standard_normal(nglobal)
    coefficients = np.ascontiguousarray(global_w[dofs])
    output = np.zeros(nglobal)
    elapsed = assemble_linear_gpu(
        module,
        layout,
        "assemble_linear",
        coords,
        coefficients,
        dofs,
        output,
        backend=backend,
        chip=chip,
    )
    if n == 1:
        reference = np.zeros_like(output)
        for xyz, local_w, indices in zip(coords, coefficients, dofs):
            np.add.at(reference, indices, _reference_stiffness_cell(xyz, degree) @ local_w)
        np.testing.assert_allclose(output, reference, rtol=1e-9, atol=1e-8)
    # grad(sum_i v_i)=0, hence the assembled vector sums to zero.
    assert abs(output.sum()) < 1e-10 * (1 + np.linalg.norm(output, 1))
    print(
        f"P{degree} {backend}: n={n}, {len(cells)} cells, {nglobal} global DOFs, "
        f"block={layout.block_shape}, grid={layout.grid_shape(len(cells))}, "
        f"quadrature points/cell={layout.quadrature_points}\n"
        f"  launch + sync: {elapsed * 1e3:.6f} ms, {nglobal / elapsed:.3e} DOFs/s"
    )


def main() -> None:
    """Check a tiny mesh, then run the selected problem size."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, choices=[1, 2], default=1)
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--backend", choices=["cuda", "amd"], default="cuda")
    parser.add_argument("--chip")
    parser.add_argument("--cells-per-block", type=int)
    args = parser.parse_args()
    if args.n < 1:
        parser.error("--n must be positive")
    chip = args.chip or ("sm_80" if args.backend == "cuda" else "gfx1100")
    assemble(1, args.degree, args.backend, chip, args.cells_per_block)
    if args.n != 1:
        assemble(args.n, args.degree, args.backend, chip, args.cells_per_block)
    print("Reference and partition-of-unity checks passed. Timing excludes setup and transfers.")


if __name__ == "__main__":
    main()
