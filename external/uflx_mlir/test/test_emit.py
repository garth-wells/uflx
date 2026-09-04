"""Correctness test for uflx_mlir.emit.generate_mlir_module.

Requires real `mlir` Python bindings (built from LLVM/MLIR source with
`-DMLIR_ENABLE_BINDINGS_PYTHON=ON` -- see README.md) importable on
PYTHONPATH; skipped entirely otherwise, since CI does not build LLVM from
source for this experimental backend.
"""

from __future__ import annotations

import ctypes

import basix
import numpy as np
import pytest
from basix_uflx import element
from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner

pytest.importorskip("mlir.ir")

from mlir.execution_engine import ExecutionEngine
from mlir.ir import Module
from mlir.passmanager import PassManager
from mlir.runtime import get_ranked_memref_descriptor
from uflx.expressions import RealScalar
from uflx_codegeneration.nodes import ArrayEntry

from uflx_mlir import generate_mlir_module, geometry_kernel_name
from uflx_mlir.emit import _alpha_signature

# Lowers scf.for/math/arith/memref/func all the way to the llvm dialect --
# the same pipeline this generator's output was validated against on real
# hardware (see the companion mlir-kernels prototype this backend
# originated from).
_PIPELINE = (
    "builtin.module("
    "convert-scf-to-cf,"
    "convert-math-to-llvm,"
    "convert-arith-to-llvm,"
    "finalize-memref-to-llvm,"
    "convert-cf-to-llvm,"
    "convert-func-to-llvm,"
    "reconcile-unrealized-casts"
    ")"
)


def _build_stiffness_form(degree: int):
    """Build `inner(grad(u), grad(v)) * dx` on a P{degree} Lagrange tetrahedron space.

    Args:
        degree: The Lagrange degree of the trial/test space.

    Returns:
        A tuple (form, ndofs).
    """
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    u = TrialFunction(space)
    v = TestFunction(space)
    return inner(grad(u), grad(v)) * dx, e.dim


def _reference_stiffness(coords: np.ndarray, degree: int) -> np.ndarray:
    """Independent, basix-quadrature-based numpy reference for the same form.

    Args:
        coords: 4x3 array of tetrahedron vertex coordinates.
        degree: The Lagrange degree of the space.

    Returns:
        The ndofs x ndofs local stiffness matrix.
    """
    e = basix.create_element(
        basix.ElementFamily.P,
        basix.CellType.tetrahedron,
        degree,
        basix.LagrangeVariant.equispaced,
    )
    qdeg = max(2 * (degree - 1), 1)
    points, weights = basix.make_quadrature(basix.CellType.tetrahedron, qdeg)
    points = np.asarray(points, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    tab = np.asarray(e.tabulate(1, points))

    x0, x1, x2, x3 = coords
    jacobian = np.column_stack([x1 - x0, x2 - x0, x3 - x0])
    det_j = np.linalg.det(jacobian)
    jacobian_inv = np.linalg.inv(jacobian)

    ndofs = tab.shape[2]
    local_matrix = np.zeros((ndofs, ndofs))
    for q, w in enumerate(weights):
        grads_ref = tab[1:4, q, :, 0].T
        grads_phys = grads_ref @ jacobian_inv
        local_matrix += (grads_phys @ grads_phys.T) * abs(det_j) * w
    return local_matrix


def _reference_geometry(coords: np.ndarray) -> np.ndarray:
    """Compute the packed affine tetrahedral Poisson metric."""
    x0, x1, x2, x3 = coords
    jacobian = np.column_stack([x1 - x0, x2 - x0, x3 - x0])
    jacobian_inv = np.linalg.inv(jacobian)
    metric = abs(np.linalg.det(jacobian)) * jacobian_inv @ jacobian_inv.T
    return metric[np.triu_indices(3)]


def _call_kernel(module: Module, kernel_name: str, a: np.ndarray, coords: np.ndarray) -> None:
    """JIT `module` and invoke `kernel_name` via MLIR's packed calling convention.

    Args:
        module: An `mlir.ir.Module` built by `generate_mlir_module`.
        kernel_name: The symbol name of the kernel function inside it.
        a: The local tensor output buffer, written in place (the kernel
            zero-initializes it itself).
        coords: The cell's coordinate dofs, read by the kernel.
    """
    with module.context:
        pm = PassManager.parse(_PIPELINE)
        pm.run(module.operation)
        engine = ExecutionEngine(module, opt_level=3)

    raw_fn = engine.lookup(kernel_name)
    a_desc_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(a)))
    coords_desc_pp = ctypes.pointer(ctypes.pointer(get_ranked_memref_descriptor(coords)))
    packed = (ctypes.c_void_p * 2)(
        ctypes.cast(a_desc_pp, ctypes.c_void_p).value,
        ctypes.cast(coords_desc_pp, ctypes.c_void_p).value,
    )
    raw_fn(packed)


def test_alpha_signature_matches_equal_dof_expressions() -> None:
    """Equivalent test/trial table expressions should share fission scratch."""
    test_expr = ArrayEntry("FE", (1, "q", "test", 0)) * RealScalar(2.0)
    trial_expr = ArrayEntry("FE", (1, "q", "trial", 0)) * RealScalar(2.0)

    assert _alpha_signature(test_expr, {"test": "$dof"}) == _alpha_signature(
        trial_expr, {"trial": "$dof"}
    )
    assert _alpha_signature(test_expr, {"test": "$dof"}) != _alpha_signature(
        trial_expr, {"q": "$dof"}
    )


@pytest.mark.parametrize("degree", [1, 2])
def test_stiffness_matches_quadrature_reference(degree: int) -> None:
    """Check generate_mlir_module's op-builder output against an independent reference.

    Builds the kernel, JITs it, calls it, and compares against a
    basix-quadrature-based numpy computation of the same form.
    """
    kernel_name = f"tabulate_tensor_p{degree}_stiffness_test"
    form, ndofs = _build_stiffness_form(degree)
    module = generate_mlir_module(
        form, degree=degree, kernel_name=kernel_name, cell=basix.CellType.tetrahedron
    )
    module_text = str(module)
    assert f"func.func @{geometry_kernel_name(kernel_name)}" in module_text
    assert f"call @{geometry_kernel_name(kernel_name)}" in module_text
    assert "memref.alloc" not in module_text.replace("memref.alloca", "")
    assert "memref<4x3xf64>" in module_text
    assert "memref<12x3xf64>" not in module_text

    coords = np.array(
        [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
        dtype=np.float64,
    )
    a = np.zeros((ndofs, ndofs), dtype=np.float64)
    _call_kernel(module, kernel_name, a, coords)

    a_ref = _reference_stiffness(coords, degree)
    np.testing.assert_allclose(a, a_ref, rtol=1e-9, atol=1e-8)

    # A constant field has zero gradient, so each row of a stiffness matrix
    # must sum to ~0 -- a cheap sanity check independent of the reference
    # computation above.
    np.testing.assert_allclose(a.sum(axis=1), 0.0, atol=1e-8)


def test_geometry_kernel_matches_numpy_reference() -> None:
    """Check the separately callable cellwise metric kernel."""
    kernel_name = "tabulate_tensor_geometry_test"
    form, _ = _build_stiffness_form(1)
    module = generate_mlir_module(
        form, degree=1, kernel_name=kernel_name, cell=basix.CellType.tetrahedron
    )

    coords = np.array(
        [[0.0, 0.3, 0.1], [1.1, -0.1, 0.05], [0.2, 1.0, -0.05], [0.15, 0.05, 1.05]],
        dtype=np.float64,
    )
    geometry = np.zeros(6, dtype=np.float64)
    _call_kernel(module, geometry_kernel_name(kernel_name), geometry, coords)

    np.testing.assert_allclose(geometry, _reference_geometry(coords), rtol=1e-12, atol=1e-12)
