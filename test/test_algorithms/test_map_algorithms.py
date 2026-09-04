"""Test map algorithms."""

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner
from uflx.functions import AbstractFunction, AbstractPhysicalFunction, AbstractReferenceFunction
from uflx.graphs.algorithms import pull_back_to_reference


def test_mass_matrix(lagrange_element):
    """Test a mass matrix."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    pulled_form = pull_back_to_reference(form.graph).root

    functions = [node for node in form.graph if isinstance(node, AbstractFunction)]
    pulled_functions = [node for node in pulled_form.graph if isinstance(node, AbstractFunction)]

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

    pulled_form = pull_back_to_reference(form.graph).root

    functions = [node for node in form.graph if isinstance(node, AbstractFunction)]
    pulled_functions = [node for node in pulled_form.graph if isinstance(node, AbstractFunction)]

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

    pulled_form = pull_back_to_reference(form.graph).root

    functions = [node for node in form.graph if isinstance(node, AbstractFunction)]
    pulled_functions = [node for node in pulled_form.graph if isinstance(node, AbstractFunction)]

    assert len(functions) == 1
    assert len(pulled_functions) == 1

    assert isinstance(functions[0], AbstractPhysicalFunction)
    assert isinstance(pulled_functions[0], AbstractReferenceFunction)
