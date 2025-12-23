"""Test forms."""

from utils import LagrangeElement, triangle

from uflx import EmbeddedCell, FunctionSpace, TestFunction, TrialFunction, dx, inner


def test_simple_form():
    """Test a simple form."""
    element = LagrangeElement(triangle, 2)
    space = FunctionSpace([(EmbeddedCell(triangle), element)])
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    print(form)
