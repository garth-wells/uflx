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
from uflx.finite_elements import AbstractReferenceMappedFiniteElement


class AbstractDomain(ABC):
    """Base class for a domain."""

    @property
    @abstractmethod
    def geometric_dimension(self) -> int:
        """The dimension of the space this domain is embedded in."""

    @property
    @abstractmethod
    def toplogical_dimension(self) -> int | None:
        """The topological dimension of the domain.

        This returns None iff the domain contains entities of a mixture
        of topological dimensions.
        """


class AbstractFiniteElementDomain(AbstractDomain):
    """Base class for a domain of a finite element function."""

    @property
    @abstractmethod
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cell types in the finite element mesh."""


class AbstractCoordinateElement(AbstractFiniteElementDomain):
    """Base class for a coordinate element.

    In a coordinate element, the geometry of the domain is defined using a
    finite element.
    """

    @property
    def element(self, cell: AbstractEntity) -> AbstractReferenceMappedFiniteElement:
        """Get the element on the given cell type."""

    @property
    def is_affine_map(self) -> bool:
        """Is the reference-to-physical map of this domain affine?"""
        return all(e.cell.is_simplex and e.lagrange_superdegree == 1 for e in self.elements)


class CoordinateElement(AbstractCoordinateElement):
    """A coordinate element."""

    def __init__(self, elements: tuple[AbstractReferenceMappedFiniteElement, ...]):
        """Initialise."""
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
    def elements(self) -> tuple[AbstractReferenceMappedFiniteElement, ...]:
        """Get the elements in the domain."""
        return self._elements

    @property
    def toplogical_dimension(self) -> int | None:
        """The topological dimension of the domain.

        This returns None iff the domain contains entities of a mixture
        of topological dimensions.
        """
        dims = {c.toplogical_dimension for c in cells}
        if len(dims) == 1:
            (dim,) = dims
            return dim
        else:
            return None


def coordinate_element(
    elements: Sequence[AbstractReferenceMappedFiniteElement] | AbstractReferenceMappedFiniteElement,
):
    """Create a domain.

    Args:
        elements: The finite element(s) used to define the geometry of the cells in this domain
    """
    if isinstance(elements, AbstractReferenceMappedFiniteElement):
        elements = (elements,)
    assert len(elements[0].reference_value_shape) == 1
    (gdim,) = elements[0].reference_value_shape
    for e in elements:
        assert e.reference_value_shape == (gdim,)

    return CoordinateElement(tuple(elements))
