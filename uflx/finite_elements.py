# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Finite element.

A finite element is an object that is used to define basis functions on a single mesh entity.
The entity on which the element is defined is called the cell.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from uflx.entities import AbstractEntity
from uflx.graphs import GraphNode
from uflx.maps import AbstractReferenceMap
from uflx.points import AbstractPoint, AbstractPointInSet
from math import prod


class AbstractFiniteElement(ABC):
    """Abstract base class for a finite element.

    To make your element library compatible with UFL, you should make a
    subclass of AbstractFiniteElement and provide implementations of all
    the abstract methods and properties. All methods and properties that
    are not marked as abstract are implemented here and should not need
    to be overwritten in your subclass.
    """

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check if this element is equal to another element."""

    @property
    @abstractmethod
    def cell(self) -> AbstractEntity:
        """Return the cell that this element is defined on."""

    @abstractmethod
    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Return the shape of the value space on a physical cell."""

    def physical_value_size(self, geometric_dimension: int) -> int:
        """Return the value size of the value space on a physical cell."""
        return prod(self.physical_value_shape(geometric_dimension))

    @property
    @abstractmethod
    def lagrange_superdegree(self) -> int | None:
        """Degree of the minimum degree Lagrange space that spans this element.

        This returns the degree of the lowest degree Lagrange space such
        that the polynomial space of the Lagrange space is a superspace
        of this element's polynomial space. If this element contains
        basis functions that are not in any Lagrange space, this
        function should return None.

        Note that on a simplex cells, the polynomial space of Lagrange
        space is a complete polynomial space, but on other cells this is
        not true. For example, on quadrilateral cells, the degree 1
        Lagrange space includes the degree 2 polynomial xy.
        """

    @property
    @abstractmethod
    def dim(self) -> int:
        """The dimension of the finite element, ie the number of basis functions."""

    @abstractmethod
    def __hash__(self):
        """Hash."""

    @abstractmethod
    def __repr__(self) -> str:
        """Representation."""


class AbstractReferenceMappedFiniteElement(AbstractFiniteElement):
    """Abstract base class for a reference-mapped finite element.

    To make your element library compatible with UFL, you should make a
    subclass of AbstractFiniteElement and provide implementations of all
    the abstract methods and properties. All methods and properties that
    are not marked as abstract are implemented here and should not need
    to be overwritten in your subclass.
    """

    @property
    @abstractmethod
    def reference_value_shape(self) -> tuple[int, ...]:
        """Return the shape of the value space on the reference cell."""

    @property
    def reference_value_size(self) -> int:
        """Return the value size of the value space on the reference cell."""
        return prod(self.reference_value_shape)

    @property
    @abstractmethod
    def reference_map(self) -> AbstractReferenceMap:
        """Get the push forward and pull back map."""

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Return the shape of the value space on a physical cell."""
        return self.reference_map.physical_value_shape(geometric_dimension)
