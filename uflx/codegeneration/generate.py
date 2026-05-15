"""Code generation."""

from typing import Self

import networkx as nx
import numpy as np

from uflx.codegeneration import symbols
from uflx.codegeneration.c import GenerateC, tables_to_c
from uflx.codegeneration.nodes import AddToLocalTensor, ArrayEntry, Loop
from uflx.complex import take_real_part
from uflx.expressions import AbstractExpression
from uflx.graphs import Graph, GraphNode, RepresentedByGraph, generate_graph, is_dag
from uflx.graphs.algorithms import replace
from uflx.integrals import AbstractIntegral, AbstractMeasure, Measure, dx
from uflx.maps import PushedForward, apply_push_forwards
from uflx.quadrature import (
    QuadratureLoop,
    QuadraturePoint,
    QuadratureRule,
    QuadratureWeight,
    quadrature_rule,
)
from uflx.scalars import RealScalar
from uflx.utils import indented


class JacobianDeterminant(AbstractExpression):
    """The determinant of the Jacobian."""

    def __init__(self, domain, point):
        """Initalise."""
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self.__class__(self.domain, self.point)


def integrals_to_quadrature(
    graph: Graph,
    rules: dict[AbstractMeasure, QuadratureRule],
    variable_namer=symbols.global_variable_namer,
) -> Graph:
    """Replace integrals with quadrature."""
    from uflx.finite_elements import EvaluatedBasisFunction
    from uflx.function_spaces import AbstractReferenceMappedFunctionSpace
    from uflx.functions import Argument
    from uflx.operators import Abs, Mult

    updated_nodes: dict[GraphNode, GraphNode] = {}
    to_replace: dict[GraphNode, GraphNode] = {}

    for node in reversed(list(nx.topological_sort(graph))):
        if isinstance(node, AbstractIntegral):
            tensor_shape_components = {}

            if not isinstance(node.measure, Measure):
                raise NotImplementedError()
            if node.measure._codim != 0 or node.measure._boundary_only:
                raise NotImplementedError("Only codim 0 integrals supported for now")

            arguments = []
            for i in nx.descendants(graph, node):
                if isinstance(i, Argument):
                    arguments.append(i)
            for i in arguments:
                i_space = i.function_space
                if not isinstance(i_space, AbstractReferenceMappedFunctionSpace):
                    raise NotImplementedError(
                        "Code generation only implemented for reference mapped spaces"
                    )
                if len(i_space.elements) != 1:
                    raise NotImplementedError(
                        "Code generation currently only implemented for spaces with "
                        "exactly one element"
                    )
                tensor_shape_components[i.component] = i_space.elements[0].dim
            tensor_shape = tuple(
                tensor_shape_components.get(i, 1)
                for i in range(
                    min(tensor_shape_components.keys()), max(tensor_shape_components.keys()) + 1
                )
            )
            variables = tuple(
                "0" if component == 1 else variable_namer.variable() for component in tensor_shape
            )

            qvariable = variable_namer.variable()

            for a in arguments:
                assert isinstance(a.function_space, AbstractReferenceMappedFunctionSpace)
                point = QuadraturePoint(rules[node.measure], variables[a.component])
                to_replace[a] = PushedForward(
                    a.function_space.elements[0].map_type,
                    EvaluatedBasisFunction(a.function_space.elements[0], qvariable, point),
                )

            domain = arguments[0].function_space.domain
            for a in arguments:
                assert a.function_space.domain == domain

            integrand = Mult(
                Mult(
                    QuadratureWeight(rules[node.measure], variables[a.component]),
                    Abs(
                        JacobianDeterminant(domain, QuadraturePoint(rules[node.measure], qvariable))
                    ),
                ),
                node.integrand,
            )

            body = AddToLocalTensor(variables, tensor_shape, integrand)

            qloop = QuadratureLoop(body, rules[node.measure], qvariable)

            next: GraphNode = qloop
            for variable, count in zip(variables[::-1], tensor_shape[::-1]):
                if variable == "0":
                    continue
                assert isinstance(count, int)
                loop = Loop(variable, 0, count, next)
                next = loop

            updated_nodes[node] = next

    new_graph = generate_graph(updated_nodes.get(graph.root, graph.root))
    new_graph = replace(new_graph, to_replace)
    return new_graph


def tabulate_quadrature(
    graph,
    variable_namer=symbols.global_variable_namer,
):
    """Generate tables of values for quadrature rules."""
    table_map = {}
    tables = {}
    to_replace = {}
    for node in graph:
        if isinstance(node, QuadratureWeight):
            if node.rule not in table_map:
                name = variable_namer.quadrature_table()
                table_map[node.rule] = name
                tables[name] = np.array(node.rule.weights)
            to_replace[node] = ArrayEntry(table_map[node.rule], (node.index,))

    return tables, replace(graph, to_replace)


def expand_jacobians(
    graph,
    variable_namer=symbols.global_variable_namer,
) -> Graph:
    """Replace jacobians with evaluations of the derivatives of finite elements."""
    from uflx.domains import AbstractCoordinateElement
    from uflx.finite_elements import EvaluatedBasisFunctionDerivative
    from uflx.operators import Add, Mult, Subtract

    to_replace: dict[GraphNode, GraphNode] = {}

    for node in graph:
        if isinstance(node, JacobianDeterminant):
            if not isinstance(node.domain, AbstractCoordinateElement):
                raise NotImplementedError()
            if len(node.domain.elements) > 1:
                raise NotImplementedError()
            (element,) = node.domain.elements
            tdim = element.cell.topological_dimension
            gdim = node.domain.geometric_dimension

            if tdim == 0 and gdim == 0:
                to_replace[node] = RealScalar(1.0)
            elif tdim == 2 and gdim == 2:
                j00: AbstractExpression = Mult(
                    ArrayEntry(symbols.coordinate_dofs, (0,)),
                    EvaluatedBasisFunctionDerivative(element, 0, node.point, (1, 0)),
                )
                j01: AbstractExpression = Mult(
                    ArrayEntry(symbols.coordinate_dofs, (0,)),
                    EvaluatedBasisFunctionDerivative(element, 0, node.point, (0, 1)),
                )
                j10: AbstractExpression = Mult(
                    ArrayEntry(symbols.coordinate_dofs, (1,)),
                    EvaluatedBasisFunctionDerivative(element, 0, node.point, (1, 0)),
                )
                j11: AbstractExpression = Mult(
                    ArrayEntry(symbols.coordinate_dofs, (1,)),
                    EvaluatedBasisFunctionDerivative(element, 0, node.point, (0, 1)),
                )

                assert isinstance(element.dim, int)
                for i in range(1, element.dim):
                    j00 = Add(
                        j00,
                        Mult(
                            ArrayEntry(symbols.coordinate_dofs, (tdim * i,)),
                            EvaluatedBasisFunctionDerivative(element, i, node.point, (1, 0)),
                        ),
                    )
                    j01 = Add(
                        j01,
                        Mult(
                            ArrayEntry(symbols.coordinate_dofs, (tdim * i,)),
                            EvaluatedBasisFunctionDerivative(element, i, node.point, (0, 1)),
                        ),
                    )
                    j10 = Add(
                        j10,
                        Mult(
                            ArrayEntry(symbols.coordinate_dofs, (tdim * i + 1,)),
                            EvaluatedBasisFunctionDerivative(element, i, node.point, (1, 0)),
                        ),
                    )
                    j11 = Add(
                        j11,
                        Mult(
                            ArrayEntry(symbols.coordinate_dofs, (tdim * i + 1,)),
                            EvaluatedBasisFunctionDerivative(element, i, node.point, (0, 1)),
                        ),
                    )

                to_replace[node] = Subtract(Mult(j00, j11), Mult(j01, j10))
            else:
                raise NotImplementedError()
    return replace(graph, to_replace)


def generate(
    form: RepresentedByGraph,
    language: str = "C",
) -> tuple[str, dict[RepresentedByGraph, str]]:
    """Generate code.

    Args:
        form: The form or other object to be assembled
        language: The programming language to use

    Returns:
        Code
    """
    from uflx.finite_elements import tabulate_finite_elements

    if language != "C":
        raise NotImplementedError("Only generation of C is supported for now")

    graph = form.graph

    assert is_dag(graph)

    # TODO: get this from somewhere
    rules: dict[AbstractMeasure, QuadratureRule] = {
        dx: quadrature_rule([[1 / 6, 1 / 6], [2 / 3, 1 / 6], [1 / 6, 2 / 3]], [1 / 6, 1 / 6, 1 / 6])
    }

    graph = integrals_to_quadrature(graph, rules)
    graph = apply_push_forwards(graph)
    graph = expand_jacobians(graph)
    graph = take_real_part(graph)

    q_tables, graph = tabulate_quadrature(graph)
    fe_tables, graph = tabulate_finite_elements(graph)
    tables = {**q_tables, **fe_tables}

    code = (
        "void tabulate_tensor_f64(\n"
        f"    double* restrict {symbols.local_tensor},\n"
        f"    const double* restrict {symbols.coefficients},\n"
        f"    const double* restrict {symbols.constants},\n"
        f"    const double* restrict {symbols.coordinate_dofs},\n"
        f"    const int* restrict {symbols.entity_local_index},\n"
        f"    const uint8_t* restrict {symbols.quadrature_permutation},\n"
        f"    void* {symbols.custom_data}\n"
        ") {\n"
    )

    code += indented(tables_to_c(tables), 2)
    code += "\n\n"
    graph.print()
    assert isinstance(graph.root, GenerateC)
    code += indented(graph.root.generate_c(), 2)
    code += "\n}\n"

    signatures = {
        form: (
            "void tabulate_tensor_f64(double* restrict, const double* restrict, "
            "const double* restrict, const double* restrict, const int* restrict, "
            "const uint8_t* restrict, void*);"
        ),
    }

    return code, signatures
