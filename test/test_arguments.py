"""Test arguments."""

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, inner
from uflx.functions import Argument


def test_argument_labelling(lagrange_element):
    """Test the integral labelling of arguments."""
    element = lagrange_element("triangle", 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    form1 = inner(u, v) * dx
    form2 = inner(u, v) * dx

    assert u.integral_label is None
    assert v.integral_label is None

    assert form1.label != form2.label
    assert form1.label is not None
    assert form2.label is not None

    for node in form1.graph:
        if isinstance(node, Argument):
            assert node.integral_label == form1.label

    for node in form2.graph:
        if isinstance(node, Argument):
            assert node.integral_label == form2.label
