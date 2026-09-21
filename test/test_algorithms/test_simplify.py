"""Test forms."""

import pytest

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
from uflx.expressions import Product, MatrixProduct
from uflx.geometry import Jacobian, JacobianInverseTranspose, JacobianInverse, JacobianTranspose
from uflx.integrals import Integral
from uflx.operators import Inner


def test_add_and_subtract_integer(lagrange_element):
    """Test that adding 2 and -2 are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    expression = u + 2 - 2
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_add_and_subtract_more_integers(lagrange_element):
    """Test that adding and subtracting integers are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    expression = u + 2 + 6 - 1 - 3 - 4
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_multiply_and_divide_integer(lagrange_element):
    """Test that multiplication then division by 2 are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    expression = u * 2 / 2
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_multiply_and_divide_more_integers(lagrange_element):
    """Test that multiplication then division by integers are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    expression = u * 4 / 6 * 3 / 2
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_add_and_subtract_function(lagrange_element):
    """Test that Function and -Function are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    f = Coefficient(space)

    expression = u - f + f
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_multiply_and_divide_integer_form(lagrange_element):
    """Test that 2 and 1/2 are successfully cancelled."""
    pytest.xfail()

    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    form = inner(2 * u, v / 2) * dx
    simpler_form = simplify(form)

    assert isinstance(form, Integral)
    assert isinstance(simpler_form, Integral)

    assert isinstance(form.integrand, Product)
    assert len(form.integrand._items) > 2

    assert isinstance(simpler_form.integrand, Product)
    assert isinstance(simpler_form.integrand._items[0], TrialFunction)
    assert isinstance(simpler_form.integrand._items[1], TestFunction)


def test_multiply_and_divide_function_form(lagrange_element):
    """Test that Function and 1/Function are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    f = Coefficient(space)

    form = (u * f) * (v / f) * dx
    simpler_form = simplify(form)

    assert isinstance(form, Integral)
    assert isinstance(simpler_form, Integral)

    assert isinstance(form.integrand, Product)
    assert len(form.integrand._items) > 2

    assert isinstance(simpler_form.integrand, Product)
    if isinstance(simpler_form.integrand._items[0], TrialFunction):
        assert isinstance(simpler_form.integrand._items[1], TestFunction)
    else:
        assert isinstance(simpler_form.integrand._items[0], TestFunction)
        assert isinstance(simpler_form.integrand._items[1], TrialFunction)


@pytest.mark.parametrize("v_first", [True, False])
@pytest.mark.parametrize("inv_first", [True, False])
@pytest.mark.parametrize("transpose", [True, False])
def test_jacobian_and_inverse_matvec(lagrange_element, v_first, inv_first, transpose):
    """Test that Jacobian and inverse Jacobian are successfully cancelled."""
    element = lagrange_element("triangle", 2, (2,))
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    if transpose:
        first = JacobianTranspose(domain)
        second = JacobianInverseTranspose(domain)
    else:
        first = Jacobian(domain)
        second = JacobianInverse(domain)
    if inv_first:
        first, second = second, first
    if v_first:
        expression = v @ first @ second
    else:
        expression = first @ second @ v

    simpler_expression = simplify(expression)

    assert isinstance(expression, MatrixProduct)
    assert isinstance(simpler_expression, TestFunction)



def test_jacobian_and_inverse_form(lagrange_element):
    """Test that Jacobian and inverse Jacobian are successfully cancelled."""
    pytest.xfail()

    element = lagrange_element("triangle", 2, (2,))
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    j = Jacobian(domain)
    j_inv_t = JacobianInverseTranspose(domain)

    form = inner(j @ u, j_inv_t @ v) * dx
    simpler_form = simplify(form)

    assert isinstance(form, Integral)
    assert isinstance(simpler_form, Integral)

    assert isinstance(form.integrand, Inner)
    assert isinstance(form.integrand.first, MatrixProduct)
    assert isinstance(form.integrand.first.first, Jacobian)
    assert isinstance(form.integrand.first.second, TrialFunction)
    assert isinstance(form.integrand.second, MatrixProduct)
    assert isinstance(form.integrand.second.first, JacobianInverseTranspose)
    assert isinstance(form.integrand.second.second, TestFunction)

    assert isinstance(simpler_form.integrand, Inner)
    if isinstance(simpler_form.integrand.first, TrialFunction):
        assert isinstance(simpler_form.integrand.second, TestFunction)
    else:
        assert isinstance(simpler_form.integrand.first, TestFunction)
        assert isinstance(simpler_form.integrand.second, TrialFunction)
