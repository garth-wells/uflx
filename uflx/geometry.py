"""Geometry."""

from typing import Any, Protocol, runtime_checkable

from uflx.domains import AbstractCoordinateElement
from uflx.expressions import AbstractExpression
from uflx.finite_elements import EvaluatedReferenceBasisFunction, EvaluatedReferenceBasisFunctionDerivative
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace
from uflx.points import AbstractPoint, Point
from uflx.scalars import Integer, RealScalar


@runtime_checkable
class ExpandableGeometry(Protocol):
    """Geometry that can be expanded into components."""

    def expand_geometry(self) -> GraphNode:
        """Expand geometry."""


class SingleSpatialCoordinate(AbstractExpression):
    """A variable represneting a component of R^d."""

    def __init__(self, dimension: int, component: int):
        """Initialise."""
        self._dimension = dimension
        self._component = component

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._dimension, self._component


class SpatialCoordinate(AbstractExpression):
    """A variable on R^d."""

    def __init__(self, dimension: int):
        """Initialise."""
        self._dimension = dimension

    def __getitem__(self, component: int) -> SingleSpatialCoordinate:
        """Get item."""
        if component < 0 or component >= self._dimension:
            raise IndexError("coordinate index out of range")
        return SingleSpatialCoordinate(self._dimension, component)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._dimension,)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._dimension,)


class ReferenceToPhysical(AbstractPoint):
    """A point mapped from the reference cell to a physical cell."""

    def __init__(self, point: AbstractPoint, domain: AbstractCoordinateElement):
        """Initialise."""
        self._point = point
        self._domain = domain

    @property
    def reference_point(self) -> AbstractPoint:
        """The point on the reference."""
        return self._point

    @property
    def domain(self) -> AbstractCoordinateElement:
        """The domain."""
        return self._domain

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._point.value_shape

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        raise NotImplementedError()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._domain

    def expand_geometry(self) -> GraphNode:
        """Expand Geometry."""
        if len(self.domain.elements) != 1:
            raise NotImplementedError("Only domains with exactly on element supported for now.")
        (element,) = self.domain.elements
        (dim,) = element.reference_value_shape

        components = [Integer(0) for _ in range(dim)]
        assert isinstance(element.dim, int)
        for i in range(element.dim):
            for j, c in enumerate(components):
                components[j] += CoordinateDofComponent(i, j, dim) * EvaluatedReferenceBasisFunction(
                    element, i, self.reference_point
                )

        return Point(*components)


class PhysicalToReference(AbstractPoint):
    """A point mapped from a physical cell to the reference cell."""

    def __init__(self, point: AbstractPoint, domain: AbstractCoordinateElement):
        """Initialise."""
        self._point = point
        self._domain = domain

    @property
    def physical_point(self) -> AbstractPoint:
        """The point on the physical cell."""
        return self._point

    @property
    def domain(self) -> AbstractCoordinateElement:
        """The domain."""
        return self._domain

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._point.value_shape

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        raise NotImplementedError()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._domain


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

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> GraphNode:
        """Expand geometry."""
        if not isinstance(self.domain, AbstractCoordinateElement):
            raise NotImplementedError()
        if len(self.domain.elements) > 1:
            raise NotImplementedError()
        (element,) = self.domain.elements
        tdim = element.cell.topological_dimension
        gdim = self.domain.geometric_dimension

        if tdim == 0 and gdim == 0:
            return RealScalar(1.0)
        elif tdim == 2 and gdim == 2:
            j00: AbstractExpression = Integer(0)
            j01: AbstractExpression = Integer(0)
            j10: AbstractExpression = Integer(0)
            j11: AbstractExpression = Integer(0)

            assert isinstance(element.dim, int)
            for i in range(element.dim):
                j00 += CoordinateDofComponent(i, 0, tdim) * EvaluatedReferenceBasisFunctionDerivative(
                    element, i, self.point, (1, 0)
                )
                j01 += CoordinateDofComponent(i, 0, tdim) * EvaluatedReferenceBasisFunctionDerivative(
                    element, i, self.point, (0, 1)
                )
                j10 += CoordinateDofComponent(i, 1, tdim) * EvaluatedReferenceBasisFunctionDerivative(
                    element, i, self.point, (1, 0)
                )
                j11 += CoordinateDofComponent(i, 1, tdim) * EvaluatedReferenceBasisFunctionDerivative(
                    element, i, self.point, (0, 1)
                )

            return j00 * j11 - j01 * j10
        else:
            raise NotImplementedError()


class CoordinateDofComponent(AbstractExpression):
    """A coordinate of a coordinate DOF."""

    def __init__(self, point, component, tdim):
        """Initialise."""
        self._point = point
        self._component = component
        self._tdim = tdim

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._component, self._tdim


def expand_geometry(
    graph: Graph,
) -> Graph:
    """Replace jacobians with evaluations of the derivatives of finite elements."""
    to_replace: dict[GraphNode, GraphNode] = {}

    for node in graph:
        if isinstance(node, GraphNode) and isinstance(node, ExpandableGeometry):
            to_replace[node] = node.expand_geometry()

    return replace(graph, to_replace)
