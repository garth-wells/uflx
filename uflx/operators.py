# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from abc import abstractmethod
from typing import Self

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

    @property
    @abstractmethod
    def init_arg(self) -> GraphNode:
        """The argument used when initialising this operator."""

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.init_arg}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        arg = self.init_arg
        if arg not in replacements:
            return self
        return self.__class__(replacements.get(arg, arg))


class BinaryOperator(AbstractOperator):
    """A binary operator.

    Binary operators act on two inputs.
    """

    @property
    @abstractmethod
    def init_args(self) -> list[GraphNode]:
        """The arguments used when initialising this operator."""

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set(self.init_args)

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        first = self.init_args[0]
        second = self.init_args[1]
        if first not in replacements and second not in replacements:
            return self
        return self.__class__(replacements.get(first, first), replacements.get(second, second))


class Inner(BinaryOperator):
    """Inner product operator.

    NOTE: document what happens here with conjugates.
    """

    def __init__(self, first: GraphNode, second: GraphNode):
        """Initialise inner product operator."""
        self._first = first
        self._second = second

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def init_args(self) -> list[GraphNode]:
        """The arguments used when initialising this operator."""
        return [self._first, self._second]


def inner(a: AbstractExpression, b: AbstractExpression) -> Inner:
    """Inner product."""
    if a.value_shape != b.value_shape:
        raise ValueError("Incompatible value shapes.")

    return Inner(a, b)


class Grad(UnaryOperator):
    """Gradient operator."""

    def __init__(self, argument: GraphNode):
        """Initialise the gradient operator."""
        self._arg = argument

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._arg.function_space.domain.geometric_dimension, )

    @property
    def init_arg(self) -> GraphNode:
        """The argument used when initialising this operator."""
        return self._arg




def grad(a: AbstractExpression) -> Grad:
    """The gradient of an expression."""
    return Grad(a)
