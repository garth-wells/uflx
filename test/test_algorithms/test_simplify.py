"""Test forms."""

from uflx import (
    Coefficient,
    TestFunction,
    TrialFunction,
    coordinate_element,
    dx,
    function_space,
    inner,
)
from uflx.algorithms import simplify
from uflx.expressions import Div, MatMult, Mult
from uflx.geometry import Jacobian, JacobianInverse
from uflx.operators import Inner


def test_function_and_inverse(lagrange_element):
    """Test that Function and 1/Function are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    f = Coefficient(space)

    form = (u * f) * (v / f) * dx

    simpler_form = simplify(form)

    form.graph.print()

    assert isinstance(form.integrand, Mult)
    assert isinstance(form.integrand.first, Mult)
    assert isinstance(form.integrand.first.first, TrialFunction)
    assert isinstance(form.integrand.first.second, Coefficient)
    assert isinstance(form.integrand.second, Div)
    assert isinstance(form.integrand.second.first, TestFunction)
    assert isinstance(form.integrand.second.second, Coefficient)

    assert isinstance(simpler_form.integrand, Mult)
    assert isinstance(simpler_form.integrand.first, TrialFunction)
    assert isinstance(form.integrand.second, TestFunction)


def test_jacobian_and_inverse(lagrange_element):
    """Test that Jacobian and inverse Jacobian are successfully cancelled."""
    element = lagrange_element("triangle", 2, (2,))
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    j = Jacobian(domain)
    inverse_j = JacobianInverse(domain)

    form = inner(j @ u, inverse_j @ v) * dx

    simpler_form = simplify(form)

    assert isinstance(form.integrand, Inner)
    assert isinstance(form.integrand.first, MatMult)
    assert isinstance(form.integrand.first.first, Jacobian)
    assert isinstance(form.integrand.first.second, TrialFunction)
    assert isinstance(form.integrand.second, MatMult)
    assert isinstance(form.integrand.second.first, JacobianInverse)
    assert isinstance(form.integrand.second.second, TestFunction)

    assert isinstance(simpler_form.integrand, Inner)
    assert isinstance(simpler_form.integrand.first, TrialFunction)
    assert isinstance(form.integrand.second, TestFunction)
