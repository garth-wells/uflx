"""Push forward and pull back maps.

These maps are uses to map function values between reference cells and physical cells
"""

from abc import ABC, abstractmethod
from typing import Any, Protocol, runtime_checkable

from uflx.expressions import AbstractExpression
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace


class AbstractReferenceMap(ABC):
    """Abstract base class for reference maps."""

    @abstractmethod
    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a reference cell to a physical cell."""

    @abstractmethod
    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a physical cell to a reference cell."""

    @abstractmethod
    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Map function values from a physical cell to a reference cell."""


class IdentityReferenceMap(AbstractReferenceMap):
    """Identity map."""

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a reference cell to a physical cell."""
        return function

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a physical cell to a reference cell."""
        return function

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Map function values from a physical cell to a reference cell."""
        return ()


class BlockedReferenceMap(AbstractReferenceMap):
    """Map for blocked element."""

    def __init__(
        self,
        component_map: AbstractReferenceMap,
        shape: tuple[int, ...],
    ):
        """Initialise."""
        self._component_map = component_map
        self._shape = shape

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a reference cell to a physical cell."""
        return function  # TODO

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a physical cell to a reference cell."""
        return function  # TODO

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Map function values from a physical cell to a reference cell."""
        return self._shape


class SymmetricReferenceMap(AbstractReferenceMap):
    """Symmetric map."""

    def __init__(
        self,
        component_map: AbstractReferenceMap,
        shape: tuple[int, ...],
        symmetry_map: dict[tuple[int, ...], int],
    ):
        """Initialise."""
        self._component_map = component_map
        self._shape = shape
        self._symmetry_map = symmetry_map

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a reference cell to a physical cell."""
        return function  # TODO

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a physical cell to a reference cell."""
        return function  # TODO

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Map function values from a physical cell to a reference cell."""
        return self._shape


class MixedReferenceMap(AbstractReferenceMap):
    """Map for a mixed element."""

    def __init__(
        self,
        sub_maps: list[AbstractReferenceMap],
        shapes: list[tuple[int, ...]],
    ):
        """Initialise."""
        self._sub_maps = sub_maps
        self._shapes = shapes

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a reference cell to a physical cell."""
        return function  # TODO

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map function values from a physical cell to a reference cell."""
        return function  # TODO

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Map function values from a physical cell to a reference cell."""
        shape: tuple[int, ...] = ()
        for s in self._shapes:
            shape += s
        return shape


@runtime_checkable
class IsPushedForward(Protocol):
    """An object that has been pushed forward."""

    def apply_push_forward(self) -> GraphNode:
        """Apply the push forward."""


@runtime_checkable
class IsPulledBack(Protocol):
    """An object that has been pulled back."""

    def apply_pull_back(self) -> GraphNode:
        """Apply the pull back."""


class PushedForward(AbstractExpression):
    """A function on a reference cell that has been mapped to a physical cell."""

    def __init__(self, map: AbstractReferenceMap, function: AbstractExpression):
        """Initialise."""
        self.map = map
        self.function = function

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.function.value_shape

    def __repr__(self):
        """Representation."""
        return f"PushedForward({self.map})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.function}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.map, self.function

    def apply_push_forward(self) -> AbstractExpression:
        """Apply the push forward."""
        return self.map.push_forward(self.function)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError("Cannot take a 'component' of a PushedForward object")


class PulledBack(AbstractExpression):
    """A function on a physical cell that has been mapped to a reference cell."""

    def __init__(self, map: AbstractReferenceMap, function: AbstractExpression):
        """Initalise."""
        self.map = map
        self.function = function

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.function.value_shape

    def __repr__(self):
        """Representation."""
        return f"PulledBack({self.map})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.function}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.map, self.function

    def apply_pull_back(self) -> AbstractExpression:
        """Apply the pull back."""
        return self.map.pull_back(self.function)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError("Cannot take a 'component' of a PulledBack object")


def apply_push_forwards(
    graph: Graph,
) -> Graph:
    """Apply push forward maps to functions."""
    return replace(
        graph,
        {
            node: node.apply_push_forward()
            for node in graph
            if isinstance(node, IsPushedForward) and isinstance(node, GraphNode)
        },
    )


def apply_pull_backs(
    graph: Graph,
) -> Graph:
    """Apply pull back maps to functions."""
    return replace(
        graph,
        {
            node: node.apply_pull_back()
            for node in graph
            if isinstance(node, IsPulledBack) and isinstance(node, GraphNode)
        },
    )
