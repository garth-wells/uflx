"""Test code generation for DG0 (cellwise-constant) Coefficients.

A DG0 Coefficient never needs a runtime dof-summation for its gradient --
grad() of a cellwise-constant physical function is resolved to an exact
zero before it ever reaches code generation (see
AbstractFunction.is_cellwise_constant and Grad.pull_back_to_reference in
uflx/operators.py), so these tests check things a purely algebraic test
(test/test_dg0_coefficients.py) cannot: that the *runtime* value of a DG0
coefficient, once summed from its one dof and carried through an actual
quadrature loop, is correct, and that a DG0 coefficient referenced only
through grad() disappears entirely from the generated ``coefficients``
array -- it costs no dofs and no summation loop, rather than merely being
present but multiplied by zero.
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


@pytest.mark.parametrize("test_degree", [1, 2])
def test_dg0_coefficient_value_vector(lagrange_element, code_dir, test_degree):
    """inner(c, v)*dx for a DG0 c equals the P0-by-Pk matrix applied to c's one dof.

    c has exactly one degree of freedom (P0), so this is the cross-space
    analogue of test_coefficient_assembly.test_coefficient_mass_vector: the
    reference matrix comes from inner(TrialFunction(P0), v)*dx, which has
    shape (v's ndofs, 1), and inner(c, v)*dx must equal that matrix applied
    to c's single dof.
    """
    p0 = lagrange_element("triangle", 0)
    v_element = lagrange_element("triangle", test_degree)
    geometry = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space0 = function_space(geometry, p0)
    v_space = function_space(geometry, v_element)
    u0 = TrialFunction(space0)
    v = TestFunction(v_space)
    c = Coefficient(space0)
    assert c.is_cellwise_constant

    ffi_ref, lib_ref = _compile(
        inner(u0, v) * dx, f"test_dg0_value_ref_degree{test_degree}", code_dir
    )
    ffi_c, lib_c = _compile(inner(c, v) * dx, f"test_dg0_value_degree{test_degree}", code_dir)

    rng = np.random.default_rng(0)
    coords = rng.random((3, 2))
    ndofs = v_element.dim

    mat = _tabulate(ffi_ref, lib_ref, (ndofs, 1), coords, None)
    for _ in range(3):
        c_dof = rng.random(1)
        vec = _tabulate(ffi_c, lib_c, (ndofs,), coords, c_dof)
        assert np.allclose(vec, mat @ c_dof)


def test_dg0_coefficient_gradient_is_exactly_zero(lagrange_element, code_dir):
    """inner(grad(c), grad(v))*dx for a DG0 c is exactly zero at runtime.

    This is the codegen-level counterpart of
    test/test_dg0_coefficients.py::test_dg0_value_and_gradient, which only
    checks that no Coefficient node survives pull_back_to_reference. Here
    the generated C is actually compiled and run, with random coordinates,
    to confirm the runtime result -- not just the graph shape -- is zero.
    """
    p0 = lagrange_element("triangle", 0)
    v_element = lagrange_element("triangle", 1)
    geometry = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space0 = function_space(geometry, p0)
    v_space = function_space(geometry, v_element)
    v = TestFunction(v_space)
    c = Coefficient(space0)

    ffi, lib = _compile(inner(grad(c), grad(v)) * dx, "test_dg0_gradient_is_zero", code_dir)

    rng = np.random.default_rng(1)
    for _ in range(3):
        coords = rng.random((3, 2))
        vec = _tabulate(ffi, lib, (v_element.dim,), coords, None)
        assert np.allclose(vec, 0.0)


def test_dg0_coefficient_absent_from_gradient_only_offsets(lagrange_element, code_dir):
    """A DG0 coefficient used only through grad() consumes no dofs at all.

    Because Grad.pull_back_to_reference (and the plain grad() function)
    resolve a cellwise-constant physical function's gradient to a zero
    before any EvaluatedReferenceCoefficientBasisFunction is ever created
    for it (see uflx/operators.py), such a coefficient never reaches
    insert_coefficient_functions's dof-offset allocation
    (external/codegeneration/uflx_codegeneration/algorithms/coefficients.py)
    -- so adding grad(c) to a form alongside a genuine P1 coefficient w1
    must not shift w1's offset or otherwise change the result, even though
    c never appears in the ``coefficients`` array passed at runtime.
    """
    p0 = lagrange_element("triangle", 0)
    p1 = lagrange_element("triangle", 1)
    geometry = coordinate_element(lagrange_element("triangle", 1, (2,)))
    v_space = function_space(geometry, p1)
    v = TestFunction(v_space)
    w1 = Coefficient(v_space)
    c = Coefficient(function_space(geometry, p0))
    assert w1.label != c.label

    ffi_plain, lib_plain = _compile(
        inner(grad(w1), grad(v)) * dx, "test_dg0_offsets_plain", code_dir
    )
    ffi_combined, lib_combined = _compile(
        inner(grad(w1) + grad(c), grad(v)) * dx, "test_dg0_offsets_combined", code_dir
    )

    rng = np.random.default_rng(2)
    coords = rng.random((3, 2))
    for _ in range(3):
        w1_dofs = rng.random(p1.dim)
        vec_plain = _tabulate(ffi_plain, lib_plain, (p1.dim,), coords, w1_dofs)
        vec_combined = _tabulate(ffi_combined, lib_combined, (p1.dim,), coords, w1_dofs)
        assert np.allclose(vec_plain, vec_combined)
