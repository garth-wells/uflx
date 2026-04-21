# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Functions.

A function is an item contained in a function space.
"""

from abc import abstractmethod

from uflx.expressions import AbstractExpression
from uflx.function_spaces import AbstractFunctionSpace


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


class TestFunction(AbstractFunction):
    """A test function."""

    __test__ = False

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise."""
        self._space = space

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space


class TrialFunction(AbstractFunction):
    """A trial function."""

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise."""
        self._space = space

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space
