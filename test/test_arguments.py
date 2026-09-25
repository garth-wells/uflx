"""Test arguments."""

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, inner
from uflx.functions import Argument
from uflx.integrals import Integral


def test_argument_labelling(lagrange_element):
    """Test the integral labelling of arguments."""
    element = lagrange_element("triangle", 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    form1 = inner(u, v) * dx
    form2 = inner(u, v) * dx

    assert isinstance(form1, Integral)
    assert isinstance(form2, Integral)

    assert u.variable is None
    assert v.variable is None

    assert form1.variable != form2.variable
    assert form1.variable is not None
    assert form2.variable is not None

    for node in form1.graph:
        if isinstance(node, Argument):
            assert node.variable == form1.variable

    for node in form2.graph:
        if isinstance(node, Argument):
            assert node.variable == form2.variable
