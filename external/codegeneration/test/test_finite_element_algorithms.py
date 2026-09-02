"""Test finite element algorithms."""

from uflx.basis_functions import EvaluatedReferenceBasisFunction
from uflx.expressions import expression_sum
from uflx.graphs import generate_graph

from uflx_codegeneration import symbols
from uflx_codegeneration.algorithms.finite_element import (
    _is_point_invariant,
    tabulate_finite_elements,
)
from uflx_codegeneration.nodes import ArrayEntry
from uflx_codegeneration.quadrature import QuadraturePoint, quadrature_rule
from uflx_codegeneration.utils import index as derivative_index


def _quadrature_point(variable: str = "q") -> QuadraturePoint:
    """A single quadrature point bound to a named loop variable."""
    rule = quadrature_rule([[0.25, 0.25]], [0.5])
    return QuadraturePoint(rule, variable)


def test_is_point_invariant_requires_a_derivative(lagrange_element):
    """A basis function's own value (derivative order 0) is never point-invariant."""
    element = lagrange_element("triangle", 1)
    assert not _is_point_invariant(element, (0, 0))


def test_is_point_invariant_true_for_affine_simplex(lagrange_element):
    """A degree 1 triangle element's gradient is constant across the cell."""
    element = lagrange_element("triangle", 1)
    assert _is_point_invariant(element, (1, 0))
    assert _is_point_invariant(element, (0, 1))


def test_is_point_invariant_false_for_higher_degree(lagrange_element):
    """A degree 2 triangle element's gradient still varies across the cell."""
    element = lagrange_element("triangle", 2)
    assert not _is_point_invariant(element, (1, 0))


def test_tabulate_finite_elements_hoists_only_the_derivative(lagrange_element):
    """tabulate_finite_elements gives a point-invariant node a constant index.

    Builds a value node (derivative order 0) and a gradient node (order 1) of the
    same degree 1 triangle element and the same quadrature point, and checks that
    only the gradient's ArrayEntry gets a constant point index -- the value node
    must keep varying with the quadrature loop variable, since (unlike its
    gradient) it isn't constant across the cell.
    """
    element = lagrange_element("triangle", 1)
    point = _quadrature_point("q")

    value_node = EvaluatedReferenceBasisFunction(element, 0, point)
    gradient_node = EvaluatedReferenceBasisFunction(element, 0, point, derivative=(1, 0))
    root = expression_sum([value_node, gradient_node])

    _, graph = tabulate_finite_elements(generate_graph(root), symbols.VariableNamer())

    entries = [n for n in graph.root.successors if isinstance(n, ArrayEntry)]
    assert len(entries) == 2
    by_derivative = {e.index[0]: e for e in entries}

    value_entry = by_derivative[derivative_index(0, 0)]
    gradient_entry = by_derivative[derivative_index(1, 0)]
    assert value_entry.index[1] == "q"
    assert gradient_entry.index[1] == 0
