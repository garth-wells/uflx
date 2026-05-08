# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Functions.

A function is an item contained in a function space.
"""

from abc import abstractmethod
from typing import Self

from uflx.expressions import AbstractExpression
from uflx.function_spaces import AbstractFunctionSpace
from uflx.graphs import GraphNode


class AbstractFunction(AbstractExpression):
    """Abstract base class for a function."""

    @property
    @abstractmethod
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.function_space.value_shape

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self


class Argument(AbstractFunction):
    """A function that is a dimension of the tensor to be assembled."""

    def __init__(self, space: AbstractFunctionSpace, component: int):
        """Initialise.

        Args:
            space: The function space that this function lives in
            component: The component of the finite element tensor
                       to be assembled that this function represents
        """
        self._space = space
        self._component = component

    @property
    def component(self) -> int:
        """The component of the finite element tensor that this function represents."""
        return self._component

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space


class TestFunction(Argument):
    """A test function."""

    __test__ = False

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise."""
        super().__init__(space, 0)


class TrialFunction(Argument):
    """A trial function."""

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise."""
        super().__init__(space, 1)
