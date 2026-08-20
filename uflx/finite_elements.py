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
from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode
from uflx.maps import AbstractReferenceMap
from uflx.points import AbstractPoint, PointInSet
from uflx.utils import product


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
        return product(self.physical_value_shape(geometric_dimension))

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
        return product(self.reference_value_shape)

    @property
    @abstractmethod
    def reference_map(self) -> AbstractReferenceMap:
        """Get the push forward and pull back map."""

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Return the shape of the value space on a physical cell."""
        return self.reference_map.physical_value_shape(geometric_dimension)


class AbstractEvaluatedBasisFunction(AbstractExpression):
    """Abstract base class for a basis function evaluated at a point in a set of points."""

    @property
    @abstractmethod
    def element(self) -> AbstractFiniteElement:
        """The finite element containing this basis function."""

    @property
    @abstractmethod
    def basis_index(self) -> int | str:
        """The index of the basis function."""

    @property
    @abstractmethod
    def point_index(self) -> int | str:
        """The index of the point in the set of points."""

    @property
    @abstractmethod
    def point(self) -> AbstractPoint:
        """The point at which the function is evaluated."""

    @property
    @abstractmethod
    def derivative(self) -> tuple[int, ...]:
        """The number of derivatives in each coordinate direction."""

    @property
    @abstractmethod
    def component(self) -> int:
        """The (flattened) component number of the basis function."""

    @property
    @abstractmethod
    def has_component(self) -> bool:
        """Check if this is an evaluation of a single component."""


class EvaluatedBasisFunction(AbstractEvaluatedBasisFunction):
    """A basis function evaluated at a point."""

    def __init__(
        self,
        element: AbstractFiniteElement,
        basis_index: int | str,
        point: AbstractPoint,
        derivative: tuple[int, ...] | None = None,
        component: int | None = None,
    ):
        """Initialise."""
        self._element = element
        self._basis_index = basis_index
        self._point = point
        if derivative is None:
            self._derivative = tuple(0 for _ in range(element.cell.topological_dimension))
        else:
            self._derivative = derivative
        if component is None and element.physical_value_size == 1:
            self._component: int | None = 0
        else:
            self._component = component

    @property
    def point(self) -> AbstractPoint:
        """The point at which the function is evaluated."""
        return self._point

    @property
    def element(self) -> AbstractFiniteElement:
        """The finite element containing this basis function."""
        return self._element

    @property
    def basis_index(self) -> int | str:
        """The index of the basis function."""
        return self._basis_index

    @property
    def point_index(self) -> int | str:
        """The index of the point in the set of points."""
        if isinstance(self._point, PointInSet):
            return self._point.index
        return 0

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def __repr__(self):
        """Representation."""
        repr = f"EvaluatedBasisFunction({self._element!r}, {self._basis_index}, {self._point!r}"
        if self._derivative is not None:
            repr += f", {self._derivative}"
        if self._component is not None:
            repr += f", {self._component}"
        repr += ")"
        return repr

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._element, self._basis_index, self._point

    @property
    def derivative(self) -> tuple[int, ...]:
        """The number of derivatives in each coordinate direction."""
        return self._derivative

    @property
    def component(self) -> int:
        """The (flattened) component number of the basis function."""
        if self._component is None:
            raise ValueError("EvaluatedBasisFunction is not an evaluation of a single component")
        return self._component

    @property
    def has_component(self) -> bool:
        """Check if this is an evaluation of a single component."""
        return self._component is not None
