"""Test coefficients."""

from uflx import TestFunction, coordinate_element, dx, function_space, inner
from uflx.functions import Coefficient
from uflx.integrals import Integral


def test_coefficient_labelling(lagrange_element):
    """Test the integral labelling of coefficients."""
    element = lagrange_element("triangle", 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w = Coefficient(space)
    v = TestFunction(space)

    form1 = inner(w, v) * dx
    form2 = inner(w, v) * dx

    assert isinstance(form1, Integral)
    assert isinstance(form2, Integral)

    assert w.variable is None
    assert v.variable is None

    assert form1.variable != form2.variable
    assert form1.variable is not None
    assert form2.variable is not None

    for node in form1.graph:
        if isinstance(node, Coefficient):
            assert node.variable == form1.variable
            assert node.label == w.label

    for node in form2.graph:
        if isinstance(node, Coefficient):
            assert node.variable == form2.variable
            assert node.label == w.label

    assert Coefficient(space).label != w.label


def test_coefficient_count_is_auto_generated(lagrange_element):
    """Test that distinct Coefficients on the same space get distinct counts."""
    element = lagrange_element("triangle", 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)

    w1 = Coefficient(space)
    w2 = Coefficient(space)
    assert w1.label != w2.label


def test_coefficient_function_space(lagrange_element):
    """Test that a Coefficient reports the function space it was built from."""
    element = lagrange_element("triangle", 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w = Coefficient(space)
    assert w.function_space == space
