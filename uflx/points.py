"""Sets of points."""

from abc import abstractmethod
from collections.abc import Hashable, Sequence
from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode


class AbstractPoint(AbstractExpression):
    """Base class for a single point in R^d."""

    @property
    @abstractmethod
    def dim(self) -> int:
        """The dimension of the point."""

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        (i,) = indices
        return PointComponent(self, i)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.dim,)


class Point(AbstractPoint):
    """A single point in R^d."""

    def __init__(self, components: Sequence[AbstractExpression]):
        """Initialise."""
        self._components = tuple(components)

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return len(self._components)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        (i,) = indices
        if isinstance(i, int):
            return self._components[i]
        return super().component(i)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set(self._components)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._components,)


class AbstractPointInSet(AbstractPoint):
    """Base class for a point in a set of points."""

    @property
    @abstractmethod
    def dim(self) -> int:
        """The dimension of the point."""

    @property
    @abstractmethod
    def index(self) -> int | str:
        """Get the index of the point in the set."""

    @property
    @abstractmethod
    def points_id(self) -> Hashable:
        """Get an identifier for the set of points."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.dim,)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()


class PointComponent(AbstractExpression):
    """A component of a point."""

    def __init__(self, point: AbstractPoint, component: int | str):
        """Initialise."""
        self._point = point
        self._component = component

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")

    @property
    def component_index(self) -> int | str:
        """Get the component of the point."""
        return self._component

    @property
    def point(self) -> AbstractPoint:
        """Get the point."""
        return self._point

    def __repr__(self):
        """Representation."""
        return f"PointComponent({self._point}, {self._component})"

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._component
