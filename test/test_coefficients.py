"""Test coefficients."""

from uflx import TestFunction, coordinate_element, dx, function_space, inner
from uflx.functions import Coefficient


def test_coefficient_labelling(lagrange_element):
    """Test the integral labelling of coefficients."""
    element = lagrange_element("triangle", 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w = Coefficient(space)
    v = TestFunction(space)

    form1 = inner(w, v) * dx
    form2 = inner(w, v) * dx

    assert w.integral_label is None
    assert v.integral_label is None

    assert form1.label != form2.label
    assert form1.label is not None
    assert form2.label is not None

    for node in form1.graph:
        if isinstance(node, Coefficient):
            assert node.integral_label == form1.label
            assert node.count == w.count

    for node in form2.graph:
        if isinstance(node, Coefficient):
            assert node.integral_label == form2.label
            assert node.count == w.count

    assert Coefficient(space).count == w.count + 1


def test_coefficient_count_is_auto_generated(lagrange_element):
    """Test that distinct Coefficients on the same space get distinct counts."""
    element = lagrange_element("triangle", 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)

    w1 = Coefficient(space)
    w2 = Coefficient(space)
    assert w1.count != w2.count


def test_coefficient_function_space(lagrange_element):
    """Test that a Coefficient reports the function space it was built from."""
    element = lagrange_element("triangle", 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w = Coefficient(space)
    assert w.function_space == space
