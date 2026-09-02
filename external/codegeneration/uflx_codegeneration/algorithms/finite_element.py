"""Finite element algorithms."""

from collections.abc import Hashable

import numpy as np
import numpy.typing as npt
from uflx.basis_functions import AbstractEvaluatedReferenceBasisFunction
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace

from uflx_codegeneration import symbols
from uflx_codegeneration.finite_element import AbstractFiniteElement
from uflx_codegeneration.nodes import ArrayEntry
from uflx_codegeneration.quadrature import QuadratureRule
from uflx_codegeneration.utils import index


def _is_point_invariant(element: AbstractFiniteElement, derivative: tuple[int, ...]) -> bool:
    """Whether tabulating `element` at this derivative gives the same row at every point.

    True for any derivative of order >= 1 of a degree 1 Lagrange element on a
    simplex cell (interval, triangle, tetrahedron): a degree 1 polynomial's
    gradient is a spatial constant everywhere on the cell, so its tabulated
    derivative is identical at every point -- unlike the element's own value
    (derivative order 0), which still varies linearly across the cell. This is
    deliberately narrow: on a non-simplex cell (eg a quadrilateral or
    hexahedron), even a degree 1 element's map is multilinear, not affine, so
    its derivative genuinely varies with position -- see
    AbstractFiniteElement.lagrange_superdegree's docstring for the same
    simplex-vs-tensor-product distinction. Restricting to exactly
    lagrange_superdegree == 1 (rather than generalising to any element whose
    derivative order happens to exceed its degree) keeps this a simple rule
    that's unambiguously correct rather than a speculative one.
    """
    return sum(derivative) >= 1 and element.cell.is_simplex and element.lagrange_superdegree == 1


def tabulate_finite_elements(
    graph: Graph,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, np.ndarray], Graph]:
    """Generate tables of values for finite elements that need to be evaluated."""
    table_map: dict[Hashable, str] = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    table_info: dict[str, tuple[AbstractFiniteElement, int, npt.NDArray[np.floating]]] = {}
    for node in graph:
        if isinstance(node, GraphNode) and isinstance(
            node, AbstractEvaluatedReferenceBasisFunction
        ):
            assert isinstance(node.element, AbstractFiniteElement)
            id = (node.element, node.point.points_set)
            if id in table_map:
                name = table_map[id]
            else:
                name = variable_namer.finite_element_table()
                table_map[id] = name
            if name not in table_info or sum(node.derivative) > table_info[name][1]:
                assert isinstance(node.point.points_set, QuadratureRule)
                table_info[name] = (
                    node.element,
                    sum(node.derivative),
                    node.point.points_set.points,
                )
            # A point-invariant node (see _is_point_invariant) reads the same value from
            # every row of the table's point axis, so index it with a constant instead of
            # the loop variable -- this is what lets hoisting analyses downstream (eg
            # uflx_mlir's hoist.py, which classifies a constant index as having no loop
            # dependency) recognise it as invariant under the quadrature loop, without the
            # table itself needing to change shape or any consumer needing new machinery.
            point_index: int | str = (
                0 if _is_point_invariant(node.element, node.derivative) else node.point_index
            )
            if node.component_index is None:
                to_replace[node] = ArrayEntry(
                    table_map[id], (index(*node.derivative), point_index, node.basis_index)
                )
            else:
                to_replace[node] = ArrayEntry(
                    table_map[id],
                    (
                        index(*node.derivative),
                        point_index,
                        node.basis_index,
                        node.component_index,
                    ),
                )

    tables = {
        name: element.tabulate(nderivs, points)
        for name, (element, nderivs, points) in table_info.items()
    }
    return tables, replace(graph, to_replace)
