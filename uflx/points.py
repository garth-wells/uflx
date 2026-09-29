"""Sets of points."""

from abc import abstractmethod
from collections.abc import Sequence
from typing import Any

from uflx.domains import RD, AbstractDomain
from uflx.expressions import AbstractExpression
from uflx.functions import AbstractVariable
from uflx.graphs import GraphNode


class AbstractSetOfPoints(AbstractDomain):
    """Base class for a set of points."""

    @property
    @abstractmethod
    def npoints(self) -> int:
        """The number of points in the set."""


class AbstractPoint(AbstractVariable):
    """Base class for a single point in R^d."""

    @property
    @abstractmethod
    def domain(self) -> AbstractSetOfPoints:
        """The set of points containing this point."""

    @property
    @abstractmethod
    def index(self) -> int | str:
        """The point's index in the set of points."""

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        (i,) = indices
        return PointComponent(self, i)


class Point(AbstractVariable):
    """A single point in R^d."""

    def __init__(self, components: Sequence[AbstractExpression], is_reference: bool = False):
        """Initialise."""
        self._components = tuple(components)
        self._is_reference = is_reference

    @property
    def is_reference(self) -> bool:
        """Check if this domain is on a reference cell."""
        return self._is_reference

    @property
    def domain(self) -> AbstractSetOfPoints:
        """The set of points containing this point."""
        return RD(len(self._components))

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

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, Point) and all(
            i == j for i, j in zip(self._components, other._components)
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.Point", *[hash(c) for c in self._components]))


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
