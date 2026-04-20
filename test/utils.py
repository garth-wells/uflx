# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Utilities for testing UFLx."""

from uflx.entity import AbstractEntity
from uflx.finite_element import AbstractReferenceMappedFiniteElement


class LagrangeElement(AbstractReferenceMappedFiniteElement):
    """A Lagrange element."""

    def __init__(self, cell: AbstractEntity, degree: int):
        """Create a Lagrange element."""
        self._cell = cell
        self._degree = degree

    def __eq__(self, other) -> bool:
        """Check if this element is equal to another element."""
        return (
            isinstance(
                other,
                LagrangeElement,
            )
            and self._cell == other._cell
            and self._degree == other._degree
        )

    @property
    def cell(self) -> AbstractEntity:
        """Return the cell that this element is defined on."""
        return self._cell

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        """Return the shape of the value space on the reference cell."""
        return ()

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Return the shape of the value space on the reference cell."""
        return ()

    @property
    def lagrange_superdegree(self) -> int | None:
        """Degree of the minimum degree Lagrange space that spans this element."""
        return self._degree


class Point(AbstractEntity):
    """A point."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Point)

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return 0

    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension."""
        match dim:
            case 0:
                return [self]
            case _:
                raise ValueError(f"Invalid dimension: {dim}")


class Interval(AbstractEntity):
    """An interval."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Interval)

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return 1

    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension."""
        match dim:
            case 0:
                return [Point() for _ in range(2)]
            case 1:
                return [self]
            case _:
                raise ValueError(f"Invalid dimension: {dim}")


class Triangle(AbstractEntity):
    """A triangle cell."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Triangle)

    @property
    def topological_dimension(self) -> int:
        """Topological dimension of the cell."""
        return 2

    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension."""
        match dim:
            case 0:
                return [Point() for _ in range(3)]
            case 1:
                return [Interval() for _ in range(3)]
            case 2:
                return [self]
            case _:
                raise ValueError(f"Invalid dimension: {dim}")


class Quadrilateral(AbstractEntity):
    """A quadrilateral."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Quadrilateral)

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return 2

    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension."""
        match dim:
            case 0:
                return [Point() for _ in range(4)]
            case 1:
                return [Interval() for _ in range(4)]
            case 2:
                return [self]
            case _:
                raise ValueError(f"Invalid dimension: {dim}")


class Tetrahedron(AbstractEntity):
    """A tetrahedron."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Tetrahedron)

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return 3

    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension."""
        match dim:
            case 0:
                return [Point() for _ in range(4)]
            case 1:
                return [Interval() for _ in range(6)]
            case 2:
                return [Triangle() for _ in range(4)]
            case 3:
                return [self]
            case _:
                raise ValueError(f"Invalid dimension: {dim}")


class Hexahedron(AbstractEntity):
    """A hexahedron."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Hexahedron)

    @property
    def topological_dimension(self) -> int:
        """Topological dimension of the cell."""
        return 3

    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension."""
        match dim:
            case 0:
                return [Point() for _ in range(4)]
            case 1:
                return [Interval() for _ in range(12)]
            case 2:
                return [Quadrilateral() for _ in range(6)]
            case 3:
                return [self]
            case _:
                raise ValueError(f"Invalid dimension: {dim}")


point = Point()
interval = Interval()
triangle = Triangle()
tetrahedron = Tetrahedron()
quadrilateral = Quadrilateral()
hexahedron = Hexahedron()
