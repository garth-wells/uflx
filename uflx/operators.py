# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from uflx.expressions import AbstractExpression, BinaryOperator, UnaryOperator
from uflx.functions import AbstractPhysicalFunction, AbstractReferenceFunction
from uflx.graphs import GraphNode
from uflx.tensors import Vector


class Inner(BinaryOperator):
    """Inner product operator.

    NOTE: document what happens here with conjugates.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class Grad(UnaryOperator):
    """Gradient operator."""

    def __init__(self, argument: GraphNode):
        """Initialise."""
        assert isinstance(argument, AbstractPhysicalFunction)
        self._physical_argument = argument
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._physical_argument.function_space.domain.geometric_dimension,)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError("Cannot get a 'component' of a Grad")


class ReferenceGrad(UnaryOperator):
    """Gradient operator."""

    def __init__(self, argument: GraphNode):
        """Initialise."""
        assert isinstance(argument, AbstractReferenceFunction)
        self._reference_argument = argument
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._reference_argument.domain_size,)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError(
            "Cannot get a 'component' of a ReferenceGrad. Try calling expand_geometry first"
        )

    def expand_geometry(self) -> GraphNode:
        """Expand geometry."""
        argument = self._reference_argument
        return Vector([argument.diff(i) for i in range(argument.domain_size)])


class Conj(UnaryOperator):
    """Complex conjugate operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return self.argument

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        raise NotImplementedError()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise NotImplementedError("Cannot get a 'component' of a Grad")
        return Conj(self.argument.component(*indices))


def grad(a: AbstractExpression) -> Grad:
    """The gradient of an expression."""
    return Grad(a)


def inner(a: AbstractExpression, b: AbstractExpression) -> AbstractExpression:
    """Inner product."""
    if a.value_shape != b.value_shape:
        raise ValueError("Incompatible value shapes.")

    if a.value_shape == ():
        return a * Conj(b)

    return Inner(a, b)
