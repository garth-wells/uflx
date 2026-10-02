"""Test that the full isotropic linear-elasticity weak form composes end to end.

Everything in uflx/operators.py and uflx/tensors.py added for linear elasticity
(sym/skew/tr/dev/transpose/Identity) has so far only been checked in isolation,
against literal numpy matrices (test_tensor_ops.py) or against a single detached
ReferenceGrad call. This module instead builds the real bilinear form

    a(u, v) = inner(sigma(u), sym(grad(v))) * dx

out of real vector-valued TrialFunction/TestFunction arguments and runs it
through the same algorithms the code-generation pipeline applies before
handing a form to quadrature lowering (see
external/codegeneration/uflx_codegeneration/generate.py): pull_back_to_reference,
apply_push_forwards, simplify. No codegen, no basix, and no numeric evaluation
is involved -- this is a pure-Python, shape- and graph-structure-level check
that the new tensor-op nodes (Transpose, Tr, and the sym/skew/dev compositions
built from them) survive that pipeline correctly, in particular that they are
correctly reconstructed around a pulled-back Grad via the generic
reconstruct_node mechanism, since neither Transpose nor Tr implements a custom
PullBackToReference rule of its own.
"""

import pytest

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner
from uflx.algorithms import pull_back_to_reference, simplify
from uflx.functions import Argument
from uflx.graphs import as_graph
from uflx.integrals import Integral
from uflx.maps import apply_push_forwards
from uflx.operators import Grad, Tr, Transpose, sym, tr
from uflx.tensors import Identity


def _sigma(displacement, lambda_, mu):
    """Isotropic Hooke's law, with hard-wired (literal float) Lame parameters."""
    strain = sym(grad(displacement))
    d = displacement.value_shape[0]
    return lambda_ * tr(strain) * Identity(d) + 2 * mu * strain


@pytest.mark.parametrize(("cell", "dim"), [("triangle", 2), ("tetrahedron", 3)])
def test_elasticity_bilinear_form_composes(lagrange_element, cell, dim):
    """inner(sigma(u), sym(grad(v))) * dx should build without error and stay scalar."""
    domain = coordinate_element(lagrange_element(cell, 1, (dim,)))
    space = function_space(domain, lagrange_element(cell, 1, (dim,)))

    u = TrialFunction(space)
    v = TestFunction(space)
    assert u.value_shape == (dim,)

    lambda_, mu = 1.7, 0.8
    sigma_u = _sigma(u, lambda_, mu)
    assert sigma_u.value_shape == (dim, dim)

    form = inner(sigma_u, sym(grad(v))) * dx
    assert isinstance(form, Integral)
    assert form.integrand.value_shape == ()


@pytest.mark.parametrize(("cell", "dim"), [("triangle", 2), ("tetrahedron", 3)])
def test_elasticity_bilinear_form_pulls_back_to_reference(lagrange_element, cell, dim):
    """The whole form must pull back to the reference cell and stay well-shaped.

    This is the generic-reconstruction check Gap 1 of the implementation plan
    called for: sym/tr/Identity are coordinate-free, so pulling grad(u) back to
    the reference cell should "just work" through reconstruct_node without a
    bespoke pullback rule for Transpose/Tr -- this test is the proof.
    """
    domain = coordinate_element(lagrange_element(cell, 1, (dim,)))
    space = function_space(domain, lagrange_element(cell, 1, (dim,)))

    u = TrialFunction(space)
    v = TestFunction(space)
    lambda_, mu = 1.7, 0.8

    form = inner(_sigma(u, lambda_, mu), sym(grad(v))) * dx

    pulled_back = pull_back_to_reference(form)
    assert isinstance(pulled_back, Integral)
    assert pulled_back.integrand.value_shape == ()

    pushed_forward = apply_push_forwards(pulled_back)
    simplified = simplify(pushed_forward)
    assert isinstance(simplified, Integral)
    assert simplified.integrand.value_shape == ()

    nodes = list(as_graph(simplified))

    # Grad works only on physical (non-reference) arguments; by this point in
    # the real pipeline every Grad must have been replaced by a reference-cell
    # equivalent, and every Argument must be reference-valued.
    assert not any(isinstance(n, Grad) for n in nodes)
    assert not any(isinstance(n, Argument) and not n.is_reference for n in nodes)

    # The tensor-op nodes themselves must have survived the round trip -- not
    # been silently dropped or left wrapping a stale (physical) Grad.
    assert any(isinstance(n, Tr) for n in nodes)
    assert any(isinstance(n, Transpose) for n in nodes)
