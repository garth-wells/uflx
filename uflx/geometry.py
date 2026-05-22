"""Geometry."""

from typing import Any

from uflx.domains import AbstractCoordinateElement
from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode
from uflx.points import AbstractPoint


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
