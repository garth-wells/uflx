"""Geometry algorithms."""

from typing import Any

from uflx.algorithms import replace
from uflx.geometry import Jacobian, JacobianDeterminant, JacobianInverse, expand_geometry
from uflx.graphs import GraphNode, as_graph
from uflx.tensors import Matrix

from uflx_codegeneration import symbols
from uflx_codegeneration.nodes import FunctionCall, Return, Variable


def insert_geometry_functions(
    expression: GraphNode,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, tuple[str, list[Variable], GraphNode]], GraphNode]:
    """Replace geometry nodes with calls to functions that compute geometry."""
    functions: dict[str, tuple[str, list[Variable], GraphNode]] = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    coordinate_dofs = Variable("const double* restrict", symbols.coordinate_dofs)
    for node in as_graph(expression):
        if isinstance(node, JacobianDeterminant):
            f = variable_namer.geometry_function_name()
            f_args: list[Any] = [coordinate_dofs]
            inputs = [coordinate_dofs]
            if node.point is not None and isinstance(node.point.index, str):
                inputs.append(Variable("int", node.point.index))
                f_args.append(node.point.index)
            to_replace[node] = FunctionCall(f, *f_args)
            functions[f] = ("double", inputs, Return(expand_geometry(node)))
        elif isinstance(node, Jacobian):
            f = variable_namer.geometry_function_name()
            fs = [
                [f"{f}_{i}_{j}" for j in range(node.value_shape[1])]
                for i in range(node.value_shape[0])
            ]
            inputs = [coordinate_dofs]
            f_args = [coordinate_dofs]
            if node.point is not None and isinstance(node.point.index, str):
                inputs.append(Variable("int", node.point.index))
                f_args.append(node.point.index)
            to_replace[node] = Matrix([[FunctionCall(f, *f_args) for f in row] for row in fs])
            for i, row in enumerate(fs):
                for j, f in enumerate(row):
                    functions[f] = (
                        "double",
                        inputs,
                        Return(expand_geometry(node.component(i, j))),
                    )
        elif isinstance(node, JacobianInverse):
            f = variable_namer.geometry_function_name()
            fs = [
                [f"{f}_{i}_{j}" for j in range(node.value_shape[1])]
                for i in range(node.value_shape[0])
            ]
            inputs = [coordinate_dofs]
            f_args = [coordinate_dofs]
            if node.point is not None and isinstance(node.point.index, str):
                inputs.append(Variable("int", node.point.index))
                f_args.append(node.point.index)
            to_replace[node] = Matrix([[FunctionCall(f, *f_args) for f in row] for row in fs])
            for i, row in enumerate(fs):
                for j, f in enumerate(row):
                    functions[f] = (
                        "double",
                        inputs,
                        Return(expand_geometry(node.component(i, j))),
                    )

    return functions, replace(expression, to_replace)
