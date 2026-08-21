"""Sets of points."""

from abc import abstractmethod
from collections.abc import Hashable
from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode


class AbstractPoint(AbstractExpression):
    """Base class for a single point in R^d."""

    @property
    @abstractmethod
    def dim(self) -> int:
        """The dimension of the point."""

    def get_component(self, component: int | str) -> AbstractExpression:
        """Get a component of the point."""
        return PointComponent(self, component)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.dim,)


class Point(AbstractPoint):
    """A single point in R^d."""

    def __init__(self, components: AbstractExpression):
        """Initialise."""
        self._components = components

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return len(self._components)

    def get_component(self, component: int | str) -> AbstractExpression:
        """Get a component of the point."""
        if isinstance(component, int):
            return self._components[component]
        return super().get_component(component)

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
    def dim(self) -> int:
        """The dimension of the point."""
        return self.points.shape[1]

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

    @property
    def component(self) -> int | str:
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
