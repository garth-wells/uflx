"""Test forms."""

from utils import LagrangeElement, triangle

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, inner


def test_simple_form():
    """Test a simple form."""
    element = LagrangeElement(triangle, 2)
    domain = coordinate_element(LagrangeElement(triangle, 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    print(form)
