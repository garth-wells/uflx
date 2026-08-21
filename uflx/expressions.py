# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Expression.

An expression is any algebraic expression that could be used as an integrand.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterable
from typing import Any

from uflx.graphs import GraphNode


class AbstractExpression(ABC):
    """Abstract base class for expressions."""

    @property
    @abstractmethod
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""

    @property
    @abstractmethod
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""

    def __mul__(self, other):
        """Multiply."""
        if not isinstance(other, AbstractExpression):
            return NotImplemented
        match len(self.value_shape), len(other.value_shape):
            case (0, _):
                return Mult(self, other)
            case (_, 0):
                return Mult(other, self)
            case (2, 1):
                return MatVec(self, other)
            case _:
                return NotImplemented

    def __add__(self, other):
        """Add."""
        if isinstance(other, AbstractExpression):
            return Add(self, other)
        return NotImplemented

    def __sub__(self, other):
        """Subtract."""
        if isinstance(other, AbstractExpression):
            return Subtract(self, other)
        return NotImplemented

    def __abs__(self):
        """Absolute value."""
        return Abs(self)

    def __repr__(self) -> str:
        """Representation."""
        return self.__class__.__name__

    @abstractmethod
    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""


class UnaryOperator(AbstractExpression):
    """A unary operator.

    Unary operators act on a single input.
    """

    def __init__(self, argument: AbstractExpression):
        """Initialise."""
        self.argument = argument

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.argument}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.argument,)

    def __repr__(self) -> str:
        """Representation."""
        return self.__class__.__name__


class BinaryOperator(AbstractExpression):
    """A binary operator.

    Binary operators act on two inputs.
    """

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        self.first = first
        self.second = second

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.first, self.second}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.first, self.second

    def __repr__(self) -> str:
        """Representation."""
        return self.__class__.__name__


class Mult(BinaryOperator):
    """Scalar multiplication operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class Add(BinaryOperator):
    """Addition operator."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        assert first.value_shape == second.value_shape
        super().__init__(first, second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.first.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return self.first.component(*indices) + self.second.component(*indices)


class Subtract(BinaryOperator):
    """Subtraction operator."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        assert first.value_shape == second.value_shape
        super().__init__(first, second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.first.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return self.first.component(*indices) - self.second.component(*indices)


class Abs(UnaryOperator):
    """Absolute value operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return Abs(self.argument.component(*indices))


class MatVec(BinaryOperator):
    """Matrix-vector multiplication operator."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        assert (
            len(first.value_shape) == 2
            and len(second.value_shape) == 1
            and first.value_shape[0] == second.value_shape[0]
        )
        super().__init__(first, second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.second.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        (index,) = indices
        return expression_sum(
            self.first.component(index, i) * self.second.component(i)
            for i in range(self.first.value_shape[1])
        )


def expression_sum(
    expressions: Iterable[AbstractExpression], default: AbstractExpression | None = None
):
    """Take the sum of a sequence of expressions."""
    result = None
    for e in expressions:
        if result is None:
            result = e
        else:
            result += e
    if result is None:
        if default is None:
            raise ValueError("Cannot sum an empty sequence without a default return value")
        return default
    return result
