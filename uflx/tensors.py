"""Tensors."""

from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode


class Matrix(AbstractExpression):
    """A matrix."""

    def __init__(self, entries: list[list[GraphNode]]):
        """Initalise."""
        self._shape = (len(entries), len(entries[0]))
        for e in entries:
            assert len(e) == self._shape[1]
        self._entries = entries

    def __repr__(self):
        """Representation."""
        return f"Matrix({self._entries})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set(i for j in self._entries for i in j)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._entries,)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._shape
