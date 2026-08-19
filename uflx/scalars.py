"""Scalar values."""

from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode


class AbstractScalar(AbstractExpression):
    """Abstract base class for scalars."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()


class AbstractInteger(AbstractScalar):
    """Abstract base class for integer values."""


class RealScalar(AbstractScalar):
    """A real scalar."""

    def __init__(self, value: float):
        """Initialise."""
        self.value = value

    def __repr__(self):
        """Representation."""
        return f"{self.value}"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.value,)

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        return f"{self.value}"


class Integer(AbstractInteger):
    """An integer."""

    def __init__(self, value: int):
        """Initialise."""
        self.value = value

    def __repr__(self):
        """Representation."""
        return f"{self.value}"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.value,)

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        return f"{self.value}"
