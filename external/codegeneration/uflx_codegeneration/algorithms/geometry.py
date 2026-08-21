"""Geometry algorithms."""

from collections.abc import Hashable
from typing import Any

import numpy as np
import numpy.typing as npt
from uflx.basis_functions import AbstractEvaluatedReferenceBasisFunction
from uflx.geometry import JacobianDeterminant, expand_geometry
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace

from uflx_codegeneration import symbols
from uflx_codegeneration.finite_element import AbstractFiniteElement
from uflx_codegeneration.nodes import ArrayEntry, FunctionCall, Variable
from uflx_codegeneration.points import AbstractPointInSet
from uflx_codegeneration.utils import index


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
