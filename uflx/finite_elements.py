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
from collections.abc import Hashable
from types import MethodType
from typing import Protocol, Self, runtime_checkable

import numpy as np

from uflx.codegeneration import symbols
from uflx.codegeneration.nodes import ArrayEntry
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace
from uflx.maps import AbstractReferenceMap
from uflx.quadrature import PointInSet
from uflx.scalars import AbstractInteger


class Dimension(AbstractInteger):
    """The dimension of a finite element."""

    def __init__(self, element: AbstractFiniteElement):
        """Initialise."""
        self._e = element

    def element(self) -> AbstractFiniteElement:
        """The finite element."""
        return self._e

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self.__class__(self._e)


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
        """Return the shape of the value space on the reference cell."""

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
    def dim(self) -> int | Dimension:
        """The dimension of the finite element, ie the number of basis functions."""
        return Dimension(self)

    @abstractmethod
    def __hash__(self):
        """Hash."""


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
    @abstractmethod
    def map_type(self) -> AbstractReferenceMap:
        """Get the push forward and pull back map."""


@runtime_checkable
class TabulatableFiniteElement(Protocol):
    """A finite element that can be tabulated."""

    def tabulate(self, points: np.ndarray, derivative: tuple[int, ...]) -> np.ndarray:
        """Create table of basis function values."""


@runtime_checkable
class CanBeTabulated(Protocol):
    """A function that can be tabulated."""

    def generate_table(self) -> np.ndarray:
        """Create table of basis function values."""

    @property
    def table_id(self) -> Hashable:
        """Get the id of the table."""


def tabulate_finite_elements(
    graph: Graph,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, np.ndarray], Graph]:
    """Generate tables of values for finite elements that need to be evaluated."""
    table_map: dict[Hashable, str] = {}
    tables = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    for node in graph:
        if (
            isinstance(node, CanBeTabulated)
            and isinstance(node, GraphNode)
            and isinstance(node, AbstractEvaluatedBasisFunction)
        ):
            id = node.table_id
            if id not in table_map:
                name = variable_namer.finite_element_table()
                table_map[id] = name
                tables[name] = node.generate_table()
            to_replace[node] = ArrayEntry(table_map[id], (node.point_index, node.basis_index))

    return tables, replace(graph, to_replace)


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


class EvaluatedBasisFunction(AbstractEvaluatedBasisFunction):
    """A basis function evaluated at a point."""

    def __init__(self, element: AbstractFiniteElement, basis_index: int | str, point: PointInSet):
        """Initalise."""
        self._element = element
        self._basis_index = basis_index
        self._point = point

        if isinstance(element, TabulatableFiniteElement):

            def generate_table(self):
                return self._element.tabulate(
                    self._point.points,
                    tuple(0 for _ in range(self._element.cell.topological_dimension)),
                )

            self.generate_table = MethodType(generate_table, self)
            self.table_id = (self.element, self.point_index)

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
        return self._point.index

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def __repr__(self):
        """Representation."""
        return f"EvaluatedBasisFunction({self._element!r}, {self._basis_index}, {self._point!r})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self.__class__(self._element, self._basis_index, self._point)


class EvaluatedBasisFunctionDerivative(AbstractEvaluatedBasisFunction):
    """A derivative of a basis function evaluated at a point."""

    def __init__(
        self,
        element: AbstractFiniteElement,
        basis_index: int | str,
        point,
        derivative: tuple[int, ...],
    ):
        """Initalise."""
        self._element = element
        self._basis_index = basis_index
        self._point = point
        self._derivative = derivative

        if isinstance(element, TabulatableFiniteElement):

            def generate_table(self):
                return self._element.tabulate(self._point.points, self._derivative)

            self.generate_table = MethodType(generate_table, self)
            self.table_id = (self._element, self._point.points_id, self._derivative)

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
        return self._point.index

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def __repr__(self):
        """Representation."""
        return (
            f"EvaluatedBasisFunctionDerivative({self._element!r}, {self._basis_index}, "
            f"{self._point!r}, {self._derivative})"
        )

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self.__class__(self._element, self._basis_index, self._point, self._derivative)
