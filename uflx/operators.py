# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from abc import abstractmethod
from typing import Self

from uflx.codegeneration.c import ConvertToCCode
from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode


class AbstractOperator(AbstractExpression):
    """Abstract base class for an operator."""

    @property
    @abstractmethod
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    @abstractmethod
    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""


class UnaryOperator(AbstractOperator):
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
        if self.argument not in replacements:
            return self
        arg = replacements.get(self.argument, self.argument)
        assert isinstance(arg, AbstractExpression)
        return self.__class__(arg)


class BinaryOperator(AbstractOperator):
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
        if self.first not in replacements and self.second not in replacements:
            return self
        first = replacements.get(self.first, self.first)
        second = replacements.get(self.second, self.second)
        assert isinstance(first, AbstractExpression)
        assert isinstance(second, AbstractExpression)
        return self.__class__(first, second)


class Inner(BinaryOperator):
    """Inner product operator.

    NOTE: document what happens here with conjugates.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()


class Mult(BinaryOperator):
    """Scalar multiplication operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.first, ConvertToCCode)
        assert isinstance(self.second, ConvertToCCode)
        return self.first.generate_c(True) + " * " + self.second.generate_c(True)


class Add(BinaryOperator):
    """Addition operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.first, ConvertToCCode)
        assert isinstance(self.second, ConvertToCCode)
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
        assert isinstance(self.first, ConvertToCCode)
        assert isinstance(self.second, ConvertToCCode)
        code = self.first.generate_c() + " - " + self.second.generate_c(True)
        if bracketed:
            return f"({code})"
        else:
            return code


class Grad(UnaryOperator):
    """Gradient operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.argument.function_space.domain.geometric_dimension,)  # type: ignore


class Conj(UnaryOperator):
    """Complex conjugate operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def re(self) -> AbstractExpression:
        """Get real part."""
        return self.argument

    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        raise NotImplementedError()


class Abs(UnaryOperator):
    """Absolute value operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.argument, ConvertToCCode)
        return f"fabs({self.argument.generate_c()})"


def grad(a: AbstractExpression) -> Grad:
    """The gradient of an expression."""
    return Grad(a)


def inner(a: AbstractExpression, b: AbstractExpression) -> AbstractExpression:
    """Inner product."""
    if a.value_shape != b.value_shape:
        raise ValueError("Incompatible value shapes.")

    if a.value_shape == ():
        return Mult(a, Conj(b))

    return Inner(a, b)
