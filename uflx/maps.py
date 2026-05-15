"""Push forward and pull back maps.

These maps are uses to map function values between reference cells and physical cells
"""

from abc import ABC, abstractmethod
from typing import Protocol, Self, runtime_checkable

from uflx.expressions import AbstractExpression
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace


class AbstractReferenceMap(ABC):
    """Abstract base class for reference maps."""

    @abstractmethod
    def push_forward_symbolic(self, function):
        """Map function values from a reference cell to a physical cell."""

    @abstractmethod
    def pull_back_symbolic(self, function):
        """Map function values from a physical cell to a reference cell."""


class IdentityReferenceMap(AbstractReferenceMap):
    """Indentity map."""

    def push_forward_symbolic(self, function):
        """Map function values from a reference cell to a physical cell."""
        return function

    def pull_back_symbolic(self, function):
        """Map function values from a physical cell to a reference cell."""
        return function


@runtime_checkable
class IsPushedForward(Protocol):
    """An object that has been pushed forward."""

    def apply_push_forward(self) -> GraphNode:
        """Apply the push forward."""


class PushedForward(AbstractExpression):
    """A function on a reference cell that has been mapped to a physical cell."""

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
        return f"PushedForward({self.map})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.function}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        function = replacements.get(self.function, self.function)
        assert isinstance(function, AbstractExpression)
        return self.__class__(self.map, function)

    def apply_push_forward(self) -> AbstractExpression:
        """Apply the push forward."""
        return self.map.push_forward_symbolic(self.function)


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
