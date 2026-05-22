"""Algorithms."""

from collections.abc import Hashable
from typing import Any, Protocol, runtime_checkable

import numpy as np
from uflx.finite_elements import AbstractEvaluatedBasisFunction
from uflx.geometry import JacobianDeterminant, expand_geometry
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace

from uflx_codegeneration import symbols
from uflx_codegeneration.nodes import ArrayEntry, FunctionCall, Variable


@runtime_checkable
class CanBeTabulated(Protocol):
    """A function that can be tabulated."""

    def generate_table(self) -> np.ndarray:
        """Create table of basis function values."""

    @property
    def table_id(self) -> Hashable:
        """Get the id of the table."""


def tabulate_finite_elements(
    graph: Graph,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, np.ndarray], Graph]:
    """Generate tables of values for finite elements that need to be evaluated."""
    table_map: dict[Hashable, str] = {}
    tables = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    for node in graph:
        if (
            isinstance(node, CanBeTabulated)
            and isinstance(node, GraphNode)
            and isinstance(node, AbstractEvaluatedBasisFunction)
        ):
            id = node.table_id
            if id not in table_map:
                name = variable_namer.finite_element_table()
                table_map[id] = name
                tables[name] = node.generate_table()
            to_replace[node] = ArrayEntry(table_map[id], (node.point_index, node.basis_index))

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
