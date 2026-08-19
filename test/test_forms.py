"""Test forms."""

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, inner


def test_simple_form(lagrange_element):
    """Test a simple form."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    print(form)
    form.graph.print()
