# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Expression.

An expression is any algebraic expression that could be used as an integrand.
"""

from abc import ABC, abstractmethod
from typing import Self, Any

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
        if isinstance(other, AbstractExpression):
            return Mult(self, other)
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


class UnaryOperator(AbstractExpression):
    """A binary operator.

    Binary operators act on two inputs.
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


class Mult(BinaryOperator):
    """Scalar multiplication operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()


class Add(BinaryOperator):
    """Addition operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()


class Subtract(BinaryOperator):
    """Subtraction operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()


class Abs(UnaryOperator):
    """Absolute value operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape
