# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from abc import abstractmethod

from uflx.expressions import AbstractExpression


class AbstractOperator(AbstractExpression):
    """Abstract base class for an operator."""

    @property
    @abstractmethod
    def arguments(self) -> tuple[AbstractExpression, ...]:
        """Arguments passed to the operator."""


class Inner(AbstractOperator):
    """Inner product operator."""

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


def inner(a: AbstractExpression, b: AbstractExpression):
    """Inner product."""
    if a.value_shape != b.value_shape:
        raise ValueError("Incompatible value shapes.")

    return Inner(a, b)
