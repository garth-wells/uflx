# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from uflx.expressions import AbstractExpression, BinaryOperator, UnaryOperator


class Inner(BinaryOperator):
    """Inner product operator.

    NOTE: document what happens here with conjugates.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()


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
