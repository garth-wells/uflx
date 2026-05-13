# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Domains.

A domain is a subset of R^d over which something can be integrated.
There is no assumption that a domain only contains cells of a single type:
one could contain (eg) a mixture of triangles and quadrilaterals, or even
a mixture of (eq) tetrahedra and intervals.
"""

from abc import ABC, abstractmethod
from collections.abc import Sequence

from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractFiniteElement


class AbstractDomain(ABC):
    """Abstract base class for a domain."""

    @property
    @abstractmethod
    def geometric_dimension(self) -> int:
        """The dimension of the space this domain is embedded in."""

    @property
    @abstractmethod
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cells in the mesh."""


class Domain(AbstractDomain):
    """A set of cells embedded in R^d."""

    def __init__(self, cells: tuple[AbstractEntity, ...], gdim: int):
        """Initialise."""
        self._cells = cells
        self._gdim = gdim

    @property
    def geometric_dimension(self) -> int:
        """Dimension of the space this domain is embedded in."""
        return self._gdim

    @property
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cells in the mesh."""
        return self._cells


def domain(cells: Sequence[AbstractEntity] | AbstractEntity, gdim: int | None = None):
    """Create a domain.

    Args:
        cells: The cell or cells included in this domain
        gdim: The geometric dimension of the space in which this domain is embedded
    """
    if isinstance(cells, AbstractEntity):
        cells = (cells,)
    if gdim is None:
        gdim = max(cell.topological_dimension for cell in cells)

    return Domain(tuple(cells), gdim)


class AbstractCoordinateElement(AbstractDomain):
    @property
    @abstractmethod
    def elements(self) -> tuple[AbstractFiniteElement, ...]:
        """Get the cells in the mesh."""


class CoordinateElement(AbstractCoordinateElement):
    def __init__(self, elements: tuple[AbstractFiniteElement, ...]):
        self._elements = elements

    @property
    def geometric_dimension(self) -> int:
        """Dimension of the space this domain is embedded in."""
        return self._elements[0].reference_value_shape[0]

    @property
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cells in the domain."""
        return tuple(e.cell for e in self._elements)

    @property
    def elements(self) -> tuple[AbstractFiniteElement, ...]:
        """Get the elements in the domain."""
        return self._elements


def coordinate_element(elements: Sequence[AbstractFiniteElement] | AbstractFiniteElement):
    """Create a domain.

    Args:
        cells: The finite element(s) used to define the geometry of the cells in this domain
    """
    if isinstance(elements, AbstractFiniteElement):
        elements = (elements,)
    assert len(elements[0].reference_value_shape) == 1
    gdim, = elements[0].reference_value_shape
    for e in elements:
        assert e.reference_value_shape == (gdim, )

    return CoordinateElement(elements)
