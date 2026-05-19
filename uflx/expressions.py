# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Expression.

An expression is any algebraic expression that could be used as an integrand.
"""

from abc import ABC, abstractmethod
from typing import Self

from uflx.codegeneration.c import GenerateC
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

    @abstractmethod
    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""

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

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        arg = replacements.get(self.argument, self.argument)
        assert isinstance(arg, AbstractExpression)
        return self.__class__(arg)


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

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        first = replacements.get(self.first, self.first)
        second = replacements.get(self.second, self.second)
        assert isinstance(first, AbstractExpression)
        assert isinstance(second, AbstractExpression)
        return self.__class__(first, second)


class Mult(BinaryOperator):
    """Scalar multiplication operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.first, GenerateC)
        assert isinstance(self.second, GenerateC)
        return self.first.generate_c(True) + " * " + self.second.generate_c(True)


class Add(BinaryOperator):
    """Addition operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.first, GenerateC)
        assert isinstance(self.second, GenerateC)
        code = self.first.generate_c() + " + " + self.second.generate_c()
        if bracketed:
            return f"({code})"
        else:
            return code


class Subtract(BinaryOperator):
    """Subtraction operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.first, GenerateC)
        assert isinstance(self.second, GenerateC)
        code = self.first.generate_c() + " - " + self.second.generate_c(True)
        if bracketed:
            return f"({code})"
        else:
            return code


class Abs(UnaryOperator):
    """Absolute value operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.argument, GenerateC)
        return f"fabs({self.argument.generate_c()})"
