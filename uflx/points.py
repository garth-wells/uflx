"""Sets of points."""

from abc import abstractmethod
from collections.abc import Sequence
from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode
from uflx.utils import Infinity


class AbstractSetOfPoints:
    """Base calss for a set of points."""

    @property
    @abstractmethod
    def npoints(self) -> int | Infinity:
        """The number of points in the set."""

    @property
    @abstractmethod
    def geometric_dimension(self) -> int:
        """The dimension of each point in the set."""


class RD(AbstractSetOfPoints):
    """R^d."""

    def __init__(self, dim: int):
        """Initialise."""
        self._dim = dim

    @property
    def npoints(self) -> int | Infinity:
        """The number of points in the set."""
        return Infinity()

    @property
    def geometric_dimension(self) -> int:
        """The dimension of each point in the set."""
        return self._dim


class AbstractPoint(AbstractExpression):
    """Base class for a single point in R^d."""

    @property
    @abstractmethod
    def dim(self) -> int:
        """The dimension of the point."""

    @property
    @abstractmethod
    def points_set(self) -> AbstractSetOfPoints:
        """The set of points containing this point."""

    @property
    def index(self) -> int | str:
        """The point's index in the set of points."""
        if self.points_set.npoints is Infinity():
            raise RuntimeError("Cannot index points in an infinite set")
        raise NotImplementedError()

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
    def points_set(self) -> AbstractSetOfPoints:
        """The set of points containing this point."""
        return RD(len(self._components))

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


class IntegrationDummyPoint(AbstractPoint):
    """A point that is a dummy variable in an integral."""

    def __init__(self, dim: int):
        """Initialise."""
        self._dim = dim

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return self._dim

    @property
    def points_set(self) -> AbstractSetOfPoints:
        """The set of points containing this point."""
        raise NotImplementedError()
