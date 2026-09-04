"""Test replace algorithm."""

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space
from uflx.graphs.algorithms import replace
from uflx.integrals import Integral


def test_replace(lagrange_element):
    """Test replace algorithm."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    form = u * dx

    replaced_graph = replace(form.graph, {u: v})

    form.graph.print()
    replaced_graph.print()

    assert isinstance(form.integrand, TrialFunction)
    assert isinstance(replaced_graph.root, Integral)
    assert isinstance(replaced_graph.root.integrand, TestFunction)
