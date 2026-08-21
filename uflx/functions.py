# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Functions.

A function is an item contained in a function space.
"""
from __future__ import annotations

from abc import abstractmethod
from typing import Any

from uflx.expressions import AbstractExpression
from uflx.function_spaces import AbstractFunctionSpace
from uflx.graphs import GraphNode
from math import prod


class AbstractFunction(AbstractExpression):
    """Abstract base class for a function."""

    @property
    @abstractmethod
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""

    @abstractmethod
    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""


class AbstractPhysicalFunction(AbstractFunction):
    """Abstract base class for a function on a physical cell."""

    @property
    @abstractmethod
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.function_space.value_shape


class AbstractReferenceFunction(AbstractFunction):
    """Abstract base class for a function on a reference cell."""


class Argument(AbstractPhysicalFunction):
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
    def component_index(self) -> int:
        """The component of the finite element tensor that this function represents."""
        return self._component

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self._component

    @property
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""
        return self._space.elements[0].cell.topological_dimension

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError()


class TestFunction(Argument):
    """A test function."""

    __test__ = False

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise."""
        super().__init__(space, 0)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._space,)

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()


class TrialFunction(Argument):
    """A trial function."""

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise."""
        super().__init__(space, 1)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._space,)

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()
