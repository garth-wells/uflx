# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from abc import abstractmethod
from typing import Self

from uflx.expressions import AbstractExpression
from uflx.graphs import union, GraphNode


class AbstractOperator(AbstractExpression):
    """Abstract base class for an operator."""

    @property
    @abstractmethod
    def arguments(self) -> tuple[AbstractExpression, ...]:
        """Arguments passed to the operator."""

    @property
    @abstractmethod
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    @abstractmethod
    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""


class Inner(AbstractOperator):
    """Inner product operator.

    NOTE: document what happens here with conjugates.
    """

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise inner product operator."""
        self._first = first
        self._second = second

    @property
    def arguments(self) -> tuple[AbstractExpression, ...]:
        """Expressions passed to the operator."""
        return (self._first, self._second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._first, self._second}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        if self._first not in replacements and self._second not in replacements:
            return self
        return Inner(
            replacements.get(self._first, self._first),
            replacements.get(self._second, self._second),
        )


def inner(a: AbstractExpression, b: AbstractExpression):
    """Inner product."""
    if a.value_shape != b.value_shape:
        raise ValueError("Incompatible value shapes.")

    return Inner(a, b)
