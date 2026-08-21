"""Algorithms."""

from collections.abc import Hashable
from typing import Any

import numpy as np
import numpy.typing as npt
from uflx.finite_elements import AbstractEvaluatedReferenceBasisFunction
from uflx.geometry import JacobianDeterminant, expand_geometry
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace
from uflx.points import PointInSet

from uflx_codegeneration import symbols
from uflx_codegeneration.finite_element import AbstractFiniteElement
from uflx_codegeneration.nodes import ArrayEntry, FunctionCall, Variable
from uflx_codegeneration.utils import index


def tabulate_finite_elements(
    graph: Graph,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, np.ndarray], Graph]:
    """Generate tables of values for finite elements that need to be evaluated."""
    table_map: dict[Hashable, str] = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    table_info: dict[str, tuple[AbstractFiniteElement, int, npt.NDArray[np.floating]]] = {}
    for node in graph:
        if isinstance(node, GraphNode) and isinstance(node, AbstractEvaluatedReferenceBasisFunction):
            assert isinstance(node.point, PointInSet)
            id = (node.element, node.point.points_id)
            if id in table_map:
                name = table_map[id]
            else:
                name = variable_namer.finite_element_table()
                table_map[id] = name
            if name not in table_info or sum(node.derivative) > table_info[name][1]:
                assert isinstance(node.element, AbstractFiniteElement)
                table_info[name] = (node.element, sum(node.derivative), node.point.points)
            if node.has_component:
                to_replace[node] = ArrayEntry(
                    table_map[id],
                    (index(*node.derivative), node.point_index, node.basis_index, node.component),
                )
            else:
                to_replace[node] = ArrayEntry(
                    table_map[id], (index(*node.derivative), node.point_index, node.basis_index)
                )

    tables = {
        name: element.tabulate(nderivs, points)
        for name, (element, nderivs, points) in table_info.items()
    }
    return tables, replace(graph, to_replace)


def insert_geometry_functions(
    graph: Graph,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, tuple[str, list[Variable], Graph]], Graph]:
    """Replace geometry nodes with calls to functions that compute geometry."""
    functions = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    coordinate_dofs = Variable("const double*", symbols.coordinate_dofs)
    for node in graph:
        if isinstance(node, JacobianDeterminant):
            f = variable_namer.geometry_function_name()
            f_args: list[Any] = [coordinate_dofs]
            inputs = [coordinate_dofs]
            if isinstance(node.point.index, str):
                inputs.append(Variable("int", node.point.index))
                f_args.append(node.point.index)
            to_replace[node] = FunctionCall(f, *f_args)
            functions[f] = ("double", inputs, expand_geometry(graph.subgraph(node)))

    return functions, replace(graph, to_replace)
