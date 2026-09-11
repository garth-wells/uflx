"""Test map algorithms."""

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner
from uflx.algorithms import pull_back_to_reference
from uflx.functions import (
    AbstractFunction,
    AbstractPhysicalFunction,
    AbstractReferenceFunction,
    Coefficient,
    ReferenceCoefficient,
)
from uflx.graphs import as_graph
from uflx.integrals import Integral


def test_mass_matrix(lagrange_element):
    """Test a mass matrix."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_reference(form)
    assert isinstance(pulled_form, Integral)

    functions = [node for node in as_graph(form) if isinstance(node, AbstractFunction)]
    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]

    assert len(functions) == 2
    assert len(pulled_functions) == 2

    for f in functions:
        assert isinstance(f, AbstractPhysicalFunction)
    for f in pulled_functions:
        assert isinstance(f, AbstractReferenceFunction)


def test_stuffness_matrix(lagrange_element):
    """Test a stiffness matrix."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(grad(u), grad(v)) * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_reference(form)
    assert isinstance(pulled_form, Integral)

    functions = [node for node in as_graph(form) if isinstance(node, AbstractFunction)]
    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]

    assert len(functions) == 2
    assert len(pulled_functions) == 2

    for f in functions:
        assert isinstance(f, AbstractPhysicalFunction)
    for f in pulled_functions:
        assert isinstance(f, AbstractReferenceFunction)


def test_linear_form(lagrange_element):
    """Test a linear form."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    v = TestFunction(space)
    form = v * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_reference(form)
    assert isinstance(pulled_form, Integral)

    functions = [node for node in as_graph(form) if isinstance(node, AbstractFunction)]
    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]

    assert len(functions) == 1
    assert len(pulled_functions) == 1

    assert isinstance(functions[0], AbstractPhysicalFunction)
    assert isinstance(pulled_functions[0], AbstractReferenceFunction)


def test_coefficient_mass_matrix_like_form(lagrange_element):
    """Test that a Coefficient pulls back like an Argument does."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w = Coefficient(space)
    v = TestFunction(space)
    form = inner(w, v) * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_reference(form)
    assert isinstance(pulled_form, Integral)

    functions = [node for node in as_graph(form) if isinstance(node, AbstractFunction)]
    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]

    assert len(functions) == 2
    assert len(pulled_functions) == 2

    for f in functions:
        assert isinstance(f, AbstractPhysicalFunction)
    for f in pulled_functions:
        assert isinstance(f, AbstractReferenceFunction)

    reference_coefficients = [f for f in pulled_functions if isinstance(f, ReferenceCoefficient)]
    assert len(reference_coefficients) == 1
    assert reference_coefficients[0].count == w.count


def test_coefficient_gradient_pulls_back(lagrange_element):
    """Test that grad(Coefficient) pulls back the same way grad(Argument) does."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w = Coefficient(space)
    v = TestFunction(space)
    form = inner(grad(w), grad(v)) * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_reference(form)
    assert isinstance(pulled_form, Integral)

    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]
    reference_coefficients = [f for f in pulled_functions if isinstance(f, ReferenceCoefficient)]
    assert len(reference_coefficients) == 1
    assert reference_coefficients[0].count == w.count


def test_distinct_coefficients_stay_distinguishable_after_pull_back(lagrange_element):
    """Test that two distinct Coefficients on the same space don't collapse together."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w1 = Coefficient(space)
    w2 = Coefficient(space)
    v = TestFunction(space)
    assert w1.count != w2.count

    form = inner(w1 + w2, v) * dx
    pulled_form = pull_back_to_reference(form)

    reference_coefficients = [
        node for node in as_graph(pulled_form) if isinstance(node, ReferenceCoefficient)
    ]
    assert len(reference_coefficients) == 2
    assert {f.count for f in reference_coefficients} == {w1.count, w2.count}
