"""Tensors."""

from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode


class Matrix(AbstractExpression):
    """A matrix."""

    def __init__(self, shape: tuple[int, int], *values: GraphNode):
        """Initalise."""
        self.shape = shape
        self.values = values

    def __repr__(self):
        """Representation."""
        return f"Matrix({self.shape}, {self.values})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set(self.values)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.shape, *self.values

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.shape
