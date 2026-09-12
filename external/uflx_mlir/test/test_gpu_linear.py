"""GPU linear-form generation and numerical checks against CPU/quadrature results."""

from __future__ import annotations

import os

import basix
import numpy as np
import pytest
from basix_uflx import element
from test_emit import (
    _build_coefficient_gradient_form,
    _call_kernel,
    _call_kernel_with_coefficients,
    _reference_stiffness,
)
from uflx import Coefficient, TestFunction, coordinate_element, dx, function_space, grad, inner

from uflx_mlir.emit import generate_mlir_module
from uflx_mlir.gpu_linear import generate_linear_assembly_gpu_module
from uflx_mlir.gpu_runtime import assemble_linear_gpu

CELL = basix.CellType.tetrahedron


def _form(degree, kind):
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    v = TestFunction(space)
    if kind == "constant":
        return v * dx, e.dim, 0
    w = Coefficient(space)
    if kind == "mass":
        return inner(w, v) * dx, e.dim, e.dim
    # Two coefficients with different element sizes expose packing/offset bugs.
    zspace = function_space(
        domain, element("Lagrange", "tetrahedron", 1, lagrange_variant="equispaced")
    )
    z = Coefficient(zspace)
    return inner(grad(w) + grad(z), grad(v)) * dx, e.dim, e.dim + 4


def _mesh_data(ndofs, ncoeff, ncells=5):
    base = np.array([[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]])
    coords = np.stack(
        [base * (1.0 + 0.13 * c) + np.array([0.1 * c, -0.2 * c, 0.03 * c]) for c in range(ncells)]
    )
    coords[1, [1, 2]] = coords[1, [2, 1]]  # Opposite orientation on one cell.
    rng = np.random.default_rng(831)
    coeffs = rng.standard_normal((ncells, ncoeff))
    maps = np.stack(
        [np.arange(ndofs, dtype=np.int32) + (c % 3) * max(1, ndofs // 2) for c in range(ncells)]
    )
    return coords, coeffs, maps


def _backend():
    if os.environ.get("CI", "").lower() == "true":
        pytest.skip("requires a GPU host")
    backend = os.environ.get("UFLX_GPU_BACKEND")
    if backend not in ("cuda", "amd"):
        pytest.skip("set UFLX_GPU_BACKEND=cuda or amd to run hardware checks")
    return backend, os.environ.get("UFLX_GPU_CHIP", "sm_89" if backend == "cuda" else "gfx1100")


@pytest.mark.parametrize("degree,threads,cells", [(1, 4, 32), (2, 16, 8), (3, 32, 4), (4, 64, 2)])
def test_automatic_launch_layout(degree, threads, cells):
    """Adapt the block to the test DOF count and group low-order cells."""
    form, ndofs = _build_coefficient_gradient_form(degree)
    module, layout = generate_linear_assembly_gpu_module(form, degree, "linear_layout", CELL)
    assert layout.ndofs == ndofs
    assert layout.coefficient_size == ndofs
    assert layout.block_shape == (threads, 1, cells)
    assert layout.grid_shape(cells + 1) == (2, 1, 1)
    assert "func.call" not in str(module)
    assert "math.absf" not in str(module)
    module.operation.verify()


@pytest.mark.parametrize("degree", [1, 2, 3, 4])
@pytest.mark.parametrize("grouping", [None, 1])
def test_gradient_vector_on_gpu(degree, grouping):
    """Scatter distinct per-cell stiffness actions, including incomplete blocks."""
    backend, chip = _backend()
    form, ndofs = _build_coefficient_gradient_form(degree)
    name = "linear_gradient"
    module, layout = generate_linear_assembly_gpu_module(
        form, degree, name, CELL, cells_per_block=grouping
    )
    coords, coeffs, maps = _mesh_data(ndofs, layout.coefficient_size)
    output = np.linspace(0.1, 0.2, int(maps.max()) + 1)
    expected = output.copy()
    for xyz, w, indices in zip(coords, coeffs, maps):
        np.add.at(expected, indices, _reference_stiffness(xyz, degree) @ w)
    elapsed = assemble_linear_gpu(
        module, layout, name, coords, coeffs, maps, output, backend=backend, chip=chip
    )
    assert elapsed > 0
    np.testing.assert_allclose(output, expected, rtol=1e-9, atol=1e-8)


@pytest.mark.parametrize("kind", ["mass", "mixed", "constant"])
def test_other_linear_forms_on_gpu(kind):
    """Reuse CPU coefficient packing for value, mixed-space and coefficient-free forms."""
    backend, chip = _backend()
    form, ndofs, ncoeff = _form(2, kind)
    name = "linear_other"
    module, layout = generate_linear_assembly_gpu_module(form, 2, name, CELL)
    assert layout.coefficient_size == ncoeff
    coords, coeffs, maps = _mesh_data(ndofs, ncoeff)
    output = np.zeros(int(maps.max()) + 1)
    expected = output.copy()
    for xyz, w, indices in zip(coords, coeffs, maps):
        cpu = generate_mlir_module(form, 2, "linear_cpu", CELL)
        local = np.zeros(ndofs)
        if ncoeff:
            _call_kernel_with_coefficients(cpu, "linear_cpu", local, xyz, w)
        else:
            _call_kernel(cpu, "linear_cpu", local, xyz)
        np.add.at(expected, indices, local)
    assemble_linear_gpu(
        module, layout, name, coords, coeffs, maps, output, backend=backend, chip=chip
    )
    np.testing.assert_allclose(output, expected, rtol=1e-9, atol=1e-8)


@pytest.mark.parametrize("grouping", [0, -1, 65, 1.5])
def test_invalid_grouping(grouping):
    """Reject invalid dimensions before kernel generation or launch."""
    form, _ = _build_coefficient_gradient_form(1)
    with pytest.raises(ValueError, match="cells_per_block"):
        generate_linear_assembly_gpu_module(form, 1, "invalid", CELL, cells_per_block=grouping)


def test_invalid_inputs_and_empty_mesh():
    """Reject bad buffer layouts and avoid launching an empty mesh."""
    form, ndofs = _build_coefficient_gradient_form(1)
    module, layout = generate_linear_assembly_gpu_module(form, 1, "empty", CELL)
    coords, coeffs, maps = _mesh_data(ndofs, ndofs)
    with pytest.raises(ValueError, match="out-of-range"):
        assemble_linear_gpu(module, layout, "empty", coords, coeffs, maps, np.zeros(1))
    with pytest.raises(ValueError, match="dtype"):
        assemble_linear_gpu(
            module, layout, "empty", coords, coeffs, maps.astype(np.int64), np.zeros(10)
        )
    output = np.ones(3)
    assert (
        assemble_linear_gpu(module, layout, "empty", coords[:0], coeffs[:0], maps[:0], output) == 0
    )
    np.testing.assert_array_equal(output, np.ones(3))
