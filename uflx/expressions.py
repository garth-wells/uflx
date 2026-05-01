# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Expression.

An expression is any algebraic expression that could be used as an integrand.
"""

from abc import ABC, abstractmethod


class AbstractExpression(ABC):
    """Abstract base class for expressions."""

    @property
    @abstractmethod
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
