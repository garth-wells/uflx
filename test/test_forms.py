"""Test forms."""

from utils import LagrangeElement, triangle

from uflx import TestFunction, TrialFunction, domain, dx, function_space, inner


def test_simple_form():
    """Test a simple form."""
    element = LagrangeElement(triangle, 2)
    space = function_space(domain(triangle), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    print(form)
