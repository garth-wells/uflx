"""Scalar values."""

from typing import Any

from uflx.expressions import AbstractExpression, Add, Mult, Subtract
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

    def __add__(self, other):
        """Add."""
        if isinstance(other, AbstractExpression):
            if self.value == 0:
                return other
            return Add(self, other)
        return NotImplemented

    def __sub__(self, other):
        """Subtract."""
        if isinstance(other, AbstractExpression):
            if self.value == 0:
                return -other
            return Subtract(self, other)
        return NotImplemented

    def __mul__(self, other):
        """Multiply."""
        if isinstance(other, AbstractExpression):
            if self.value == 0:
                return self
            if self.value == 1:
                return other
            return Mult(self, other)
        return NotImplemented
