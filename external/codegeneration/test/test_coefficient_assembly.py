"""Test code generation for Coefficients.

These tests check a Coefficient's generated code against the *independently*
generated code for the corresponding bilinear form, rather than against
hand-derived reference numbers: a Coefficient's value is, by definition, a
runtime sum over its own dofs (see
``uflx_codegeneration.coefficients.EvaluatedReferenceCoefficientBasisFunction``),
so ``inner(w, v) * dx`` for a Coefficient ``w`` on the same space as ``u``
must equal the ordinary mass matrix (from ``inner(u, v) * dx``) applied to
``w``'s dof vector, and similarly for the stiffness matrix with
``inner(grad(w), grad(v)) * dx``.
"""

from typing import Any

import numpy as np
import pytest
from cffi import FFI
from uflx import (
    Coefficient,
    TestFunction,
    TrialFunction,
    coordinate_element,
    dx,
    function_space,
    grad,
    inner,
)

import uflx_codegeneration


def _compile(form, name, code_dir):
    """Generate and compile C code for a form, returning the ffi and library handle."""
    code, signature = uflx_codegeneration.generate(form)
    ffi = FFI()
    ffi.cdef(signature)
    ffi.set_source(name, code)
    so = ffi.compile(code_dir)
    return ffi, ffi.dlopen(so)


def _tabulate(ffi, lib: Any, shape: tuple[int, ...], coords: np.ndarray, w: np.ndarray | None):
    """Call tabulate_tensor_f64 and return the resulting local tensor."""
    tensor = np.zeros(shape)
    empty = np.zeros(0)
    w_array = empty if w is None else np.ascontiguousarray(w, dtype=np.float64)
    lib.tabulate_tensor_f64(
        ffi.cast("double*", tensor.ctypes.data),
        ffi.cast("double*", w_array.ctypes.data),
        ffi.cast("double*", empty.ctypes.data),
        ffi.cast("double*", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    return tensor


@pytest.mark.parametrize("degree", [1, 2])
def test_coefficient_mass_vector(lagrange_element, code_dir, degree):
    """inner(w, v) * dx should equal the mass matrix applied to w's dof vector."""
    element = lagrange_element("triangle", degree)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    w = Coefficient(space)

    ffi_mass, lib_mass = _compile(
        inner(u, v) * dx, f"test_coefficient_mass_ref_degree{degree}", code_dir
    )
    ffi_coeff, lib_coeff = _compile(
        inner(w, v) * dx, f"test_coefficient_mass_degree{degree}", code_dir
    )

    rng = np.random.default_rng(0)
    coords = rng.random((3, 2))
    ndofs = element.dim

    mass_matrix = _tabulate(ffi_mass, lib_mass, (ndofs, ndofs), coords, None)
    for _ in range(3):
        w_dofs = rng.random(ndofs)
        vec = _tabulate(ffi_coeff, lib_coeff, (ndofs,), coords, w_dofs)
        assert np.allclose(vec, mass_matrix @ w_dofs)


@pytest.mark.parametrize("degree", [1, 2])
def test_coefficient_gradient_vector(lagrange_element, code_dir, degree):
    """inner(grad(w), grad(v))*dx should equal the stiffness matrix applied to w's dofs.

    This exercises differentiating a Coefficient: the runtime dof-summation
    must be deferred until after grad() has been resolved, otherwise
    differentiating the (already-summed) value would be meaningless -- see
    EvaluatedReferenceCoefficientBasisFunction.
    """
    element = lagrange_element("triangle", degree)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    w = Coefficient(space)

    ffi_stiff, lib_stiff = _compile(
        inner(grad(u), grad(v)) * dx, f"test_coefficient_stiffness_ref_degree{degree}", code_dir
    )
    ffi_coeff, lib_coeff = _compile(
        inner(grad(w), grad(v)) * dx, f"test_coefficient_stiffness_degree{degree}", code_dir
    )

    rng = np.random.default_rng(1)
    coords = rng.random((3, 2))
    ndofs = element.dim

    stiffness_matrix = _tabulate(ffi_stiff, lib_stiff, (ndofs, ndofs), coords, None)
    for _ in range(3):
        w_dofs = rng.random(ndofs)
        vec = _tabulate(ffi_coeff, lib_coeff, (ndofs,), coords, w_dofs)
        assert np.allclose(vec, stiffness_matrix @ w_dofs)


def test_coefficient_offsets_are_distinct(lagrange_element, code_dir):
    """Two distinct Coefficients on the same space get non-overlapping offsets into w."""
    degree = 1
    element = lagrange_element("triangle", degree)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    w1 = Coefficient(space)
    w2 = Coefficient(space)
    assert w1.count != w2.count

    ffi_mass, lib_mass = _compile(inner(u, v) * dx, "test_coefficient_offset_mass_ref", code_dir)
    ffi_coeff, lib_coeff = _compile(inner(w1 + w2, v) * dx, "test_coefficient_offset_sum", code_dir)

    rng = np.random.default_rng(2)
    coords = rng.random((3, 2))
    ndofs = element.dim

    mass_matrix = _tabulate(ffi_mass, lib_mass, (ndofs, ndofs), coords, None)
    w1_dofs = rng.random(ndofs)
    w2_dofs = rng.random(ndofs)
    vec = _tabulate(ffi_coeff, lib_coeff, (ndofs,), coords, np.concatenate([w1_dofs, w2_dofs]))
    assert np.allclose(vec, mass_matrix @ (w1_dofs + w2_dofs))
