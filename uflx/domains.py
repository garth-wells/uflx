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
    """Abstract base class for a domain."""

    @property
    @abstractmethod
    def geometric_dimension(self) -> int:
        """The dimension of the space this domain is embedded in."""

    @property
    @abstractmethod
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cells in the mesh."""


class AbstractCoordinateElement(AbstractDomain):
    """Abstract coordinate element.

    In a coordinate element, the geometry of the cell is represented by a finite element.
    """

    @property
    @abstractmethod
    def elements(self) -> tuple[AbstractReferenceMappedFiniteElement, ...]:
        """Get the cells in the mesh."""

    @property
    def is_affine_map(self) -> bool:
        """Whether the reference-to-physical map of this domain is affine.

        True iff every element used to define the domain's coordinates is a degree 1
        Lagrange space on a simplex cell -- the one case where the Jacobian of the
        map is spatially constant (the same at every point of the cell) rather than
        varying with position. A degree 1 Lagrange element on a non-simplex cell (eg
        a quadrilateral or hexahedron) still gives a multilinear, not affine, map:
        see AbstractFiniteElement.lagrange_superdegree's docstring for the same
        simplex-vs-tensor-product distinction. This is a purely static property --
        it depends only on cell shape and polynomial degree, never on the mesh's
        actual coordinate values -- so callers (eg code generation deciding whether
        Jacobian-derived quantities can be hoisted out of the quadrature loop
        entirely) can rely on it without any per-cell computation.
        """
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
