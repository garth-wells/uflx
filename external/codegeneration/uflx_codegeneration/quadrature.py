"""Quadrature rules."""

from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
from uflx.algorithms import replace
from uflx.basis_functions import EvaluatedPhysicalBasisFunction, EvaluatedReferenceBasisFunction
from uflx.domains import AbstractCoordinateElement, AbstractDomain
from uflx.expressions import AbstractExpression
from uflx.function_spaces import AbstractReferenceMappedFunctionSpace
from uflx.functions import AbstractPhysicalFunction, Argument, ReferenceArgument
from uflx.geometry import (
    Jacobian,
    JacobianDeterminant,
    JacobianInverse,
    JacobianInverseTranspose,
    JacobianTranspose,
    ReferenceToPhysical,
    SingleSpatialCoordinate,
)
from uflx.graphs import Graph, GraphNode, as_graph
from uflx.integrals import AbstractIntegral, AbstractMeasure, Measure
from uflx.points import AbstractPoint, AbstractSetOfPoints, Point, PointComponent

from uflx_codegeneration import symbols
from uflx_codegeneration.c import GenerateC
from uflx_codegeneration.nodes import AddToLocalTensor, ArrayEntry, Loop
from uflx_codegeneration.utils import indented


class QuadratureRule(AbstractSetOfPoints):
    """A quadrature rule."""

    def __init__(self, points: npt.NDArray[np.floating], weights: npt.NDArray[np.floating]):
        """Initialise."""
        assert points.shape[0] == len(weights)
        self.points = points
        self.weights = weights

    @property
    def npoints(self) -> int:
        """The number of points in the set."""
        return len(self.weights)

    @property
    def geometric_dimension(self) -> int:
        """The dimension of each point in the set."""
        return self.points.shape[1]


class QuadraturePoint(AbstractPoint):
    """A point in a quadrature rule."""

    def __init__(self, rule: QuadratureRule, index: int | str):
        """Initialise."""
        self.rule = rule
        self._index = index

    @property
    def points_set(self) -> QuadratureRule:
        """Get all the points in the set."""
        return self.rule

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return self.rule.geometric_dimension

    @property
    def index(self) -> int | str:
        """Get the index of the point in the set."""
        return self._index

    def __repr__(self):
        """Representation."""
        return f"QuadraturePoint({self._index})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.rule, self._index


class QuadratureWeight(AbstractExpression):
    """A weight in a quadrature rule."""

    def __init__(self, rule: QuadratureRule, index: int | str):
        """Initialise."""
        self.rule = rule
        self._index = index

    @property
    def index(self) -> int | str:
        """Get the index of the point in the set."""
        return self._index

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def __repr__(self):
        """Representation."""
        return f"QuadratureWeight({self.index})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.rule, self._index

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class QuadratureLoop:
    """A loop over the points in a quadrature rule."""

    def __init__(self, body: GraphNode, rule: QuadratureRule, variable: str):
        """Initialise."""
        self.body = body
        self.rule = rule
        self.variable = variable

    def __repr__(self):
        """Representation."""
        return f"QuadratureLoop({self.variable})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.body, self.rule, self.variable

    def generate_c(self) -> str:
        """Generate code for this object."""
        if not isinstance(self.body, GenerateC):
            raise NotImplementedError(f"GenerateC is not implemented for {self.body.__class__}")
        return (
            f"for (int {self.variable}=0; {self.variable}!={self.rule.npoints}; "
            f"++{self.variable})\n"
            "{\n" + indented(self.body.generate_c(), 2) + "\n}"
        )


def quadrature_rule(
    points: Sequence[Sequence[float]] | npt.ArrayLike,
    weights: Sequence[float] | npt.ArrayLike,
) -> QuadratureRule:
    """Create a quadrature rule."""
    return QuadratureRule(np.array(points), np.array(weights))


def extract_domain(graph: Graph, node: GraphNode) -> AbstractDomain:
    """Extract the domain associated with a node."""
    domain: AbstractDomain | None = None
    for i in graph.descendants(node):
        if isinstance(i, AbstractPhysicalFunction):
            if domain is None:
                domain = i.function_space.domain
            else:
                assert domain == i.function_space.domain
        if hasattr(i, "domain"):
            if domain is None:
                domain = i.domain
            else:
                assert domain == i.domain
    assert domain is not None
    return domain


def integrals_to_quadrature(
    expression: GraphNode,
    rules: dict[AbstractMeasure, QuadratureRule],
    variable_namer=symbols.global_variable_namer,
) -> GraphNode:
    """Replace integrals with quadrature."""
    updated_nodes: dict[GraphNode, GraphNode] = {}
    to_replace: dict[GraphNode, GraphNode] = {}

    graph = as_graph(expression)
    for node in graph.ordered_nodes():
        if isinstance(node, AbstractIntegral):
            rule = rules[node.measure]
            qvariable = variable_namer.variable()
            qpoint = QuadraturePoint(rule, qvariable)

            tensor_shape_components = {}

            if not isinstance(node.measure, Measure):
                raise NotImplementedError()
            if node.measure._codim != 0 or node.measure._boundary_only:
                raise NotImplementedError("Only codim 0 integrals supported for now")

            arguments = []
            for i in graph.descendants(node):
                if isinstance(i, (Argument, ReferenceArgument)) and i.integral_label == node.label:
                    arguments.append(i)
                if isinstance(i, SingleSpatialCoordinate):
                    domain = extract_domain(graph, node)
                    if not isinstance(domain, AbstractCoordinateElement):
                        raise NotImplementedError(
                            "Code generation only implemented for reference mapped domain"
                        )
                    to_replace[i] = PointComponent(
                        ReferenceToPhysical(qpoint, domain), i._component
                    )
                if isinstance(i, Jacobian) and i.point is None:
                    to_replace[i] = Jacobian(i.domain, qpoint)
                if isinstance(i, JacobianInverse) and i.point is None:
                    to_replace[i] = JacobianInverse(i.domain, qpoint)
                if isinstance(i, JacobianTranspose) and i.point is None:
                    to_replace[i] = JacobianTranspose(i.domain, qpoint)
                if isinstance(i, JacobianInverseTranspose) and i.point is None:
                    to_replace[i] = JacobianInverseTranspose(i.domain, qpoint)
                if isinstance(i, JacobianDeterminant) and i.point is None:
                    to_replace[i] = JacobianDeterminant(i.domain, qpoint)
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
                tensor_shape_components[i.component_index] = i_space.elements[0].dim
            tensor_shape = tuple(
                tensor_shape_components.get(i, 1)
                for i in range(
                    min(tensor_shape_components.keys()), max(tensor_shape_components.keys()) + 1
                )
            )
            variables = tuple(
                "0" if component == 1 else variable_namer.variable() for component in tensor_shape
            )

            for a in arguments:
                assert isinstance(a.function_space, AbstractReferenceMappedFunctionSpace)
                assert isinstance(a.function_space.domain, AbstractCoordinateElement)
                if isinstance(a, Argument):
                    to_replace[a] = EvaluatedPhysicalBasisFunction(
                        a.function_space,
                        a.function_space.elements[0],
                        variables[a.component_index],
                        ReferenceToPhysical(qpoint, a.function_space.domain),
                    )
                elif isinstance(a, ReferenceArgument):
                    to_replace[a] = EvaluatedReferenceBasisFunction(
                        a.function_space.elements[0],
                        variables[a.component_index],
                        qpoint,
                    )

            domain = arguments[0].function_space.domain
            for a in arguments:
                assert a.function_space.domain == domain

            assert isinstance(domain, AbstractCoordinateElement)
            integrand = QuadratureWeight(rules[node.measure], qvariable) * node.integrand

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

    return replace(updated_nodes.get(graph.root, graph.root), to_replace)


def tabulate_quadrature(
    expression: GraphNode,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, npt.NDArray(np.floating)], GraphNode]:
    """Generate tables of values for quadrature rules."""
    table_map = {}
    tables = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    for node in as_graph(expression):
        if isinstance(node, QuadratureWeight):
            id = (node.rule, "weights")
            if id not in table_map:
                name = variable_namer.quadrature_table()
                table_map[id] = name
                tables[name] = node.rule.weights
            to_replace[node] = ArrayEntry(table_map[id], (node.index,))
        if isinstance(node, QuadraturePoint):
            id = (node.rule, "points")
            if id not in table_map:
                name = variable_namer.quadrature_table()
                table_map[id] = name
                tables[name] = node.rule.points
            to_replace[node] = Point(
                [
                    ArrayEntry(
                        table_map[id],
                        (
                            node.dim * node.index + i
                            if isinstance(node.index, int)
                            else f"{node.dim} * {node.index} + {i}",
                        ),
                    )
                    for i in range(node.dim)
                ]
            )

    return tables, replace(expression, to_replace)
