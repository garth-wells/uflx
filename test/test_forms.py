"""Test forms."""

from uflx import TestFunction, TrialFunction, FunctionSpace, EmbeddedCell, dx, inner
from utils import triangle, LagrangeElement


def test_simple_form():
    element = LagrangeElement(triangle, 2)
    space = FunctionSpace([(EmbeddedCell(triangle), element)])
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    print(form)
