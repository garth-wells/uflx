"""Test graph algorithms."""

from utils import LagrangeElement, triangle

from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space
from uflx.graphs.algorithms import replace


def test_replace():
    """Test replace algorithm."""
    element = LagrangeElement(triangle, 2)
    domain = coordinate_element(LagrangeElement(triangle, 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    form1 = u * dx
    form2 = v * dx

    replaced_graph1 = replace(form1.graph, {u: v})
    graph2 = form2.graph

    assert graph2.root == replaced_graph1.root
