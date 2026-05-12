# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Function spaces.

A function space is a space containing functions defined on a domain.
In most if not all cases, these will be finite dimensional.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence

from uflx.domains import AbstractDomain
from uflx.finite_elements import AbstractReferenceMappedFiniteElement


class AbstractFunctionSpace(ABC):
    """Abstract base class for a function space."""

    @property
    @abstractmethod
    def domain(self) -> AbstractDomain:
        """Domain of the function space."""

    @property
    @abstractmethod
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the function space."""


class AbstractReferenceMappedFunctionSpace(AbstractFunctionSpace):
    """Abstract base class for a function space whose functions are mapped from a reference cell."""

    @property
    @abstractmethod
    def elements(self) -> tuple[AbstractReferenceMappedFiniteElement, ...]:
        """Elements in the function space."""


class FunctionSpace(AbstractReferenceMappedFunctionSpace):
    """Function space."""

    def __init__(
        self,
        domain: AbstractDomain,
        elements: tuple[AbstractReferenceMappedFiniteElement, ...],
    ):
        """Initialise."""
        self._domain = domain
        self._elements = elements
        gdim = domain.geometric_dimension
        shape = elements[0].physical_value_shape(gdim)
        for element in elements[1:]:
            if element.physical_value_shape(gdim) != shape:
                raise ValueError(
                    "Elements in a functions space must have the same physical value shape."
                )

    @property
    def domain(self) -> AbstractDomain:
        """Domain of the function space."""
        return self._domain

    @property
    def elements(self) -> tuple[AbstractReferenceMappedFiniteElement, ...]:
        """Elements in the function space."""
        return self._elements

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the function space."""
        return self.elements[0].physical_value_shape(self.domain.geometric_dimension)


def function_space(
    domain: AbstractDomain,
    elements: Sequence[AbstractReferenceMappedFiniteElement] | AbstractReferenceMappedFiniteElement,
) -> FunctionSpace:
    """Create a function space.

    Args:
        domain: The domain on which the function space is defined.
        elements: The elements in the funciton space.
    """
    if isinstance(elements, AbstractReferenceMappedFiniteElement):
        elements = (elements,)
    return FunctionSpace(domain, tuple(elements))
