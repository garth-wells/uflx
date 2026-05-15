# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Utilities for testing UFLx."""

import numpy as np

from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractReferenceMappedFiniteElement, Dimension
from uflx.maps import AbstractReferenceMap, IdentityReferenceMap


class LagrangeElement(AbstractReferenceMappedFiniteElement):
    """A Lagrange element."""

    def __init__(
        self, cell: AbstractEntity, degree: int, block_shape: tuple[int, ...] | None = None
    ):
        """Create a Lagrange element."""
        self._cell = cell
        self._degree = degree
        self._block_shape = block_shape

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
        if self._block_shape is None:
            return ()
        return self._block_shape

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Return the shape of the value space on the reference cell."""
        return ()

    @property
    def lagrange_superdegree(self) -> int | None:
        """Degree of the minimum degree Lagrange space that spans this element."""
        return self._degree

    @property
    def dim(self) -> int | Dimension:
        """The dimension of the finite element, ie the number of basis functions."""
        if isinstance(self._cell, Point):
            return 1
        if isinstance(self._cell, Interval):
            return self._degree + 1
        if isinstance(self._cell, Triangle):
            return (self._degree + 1) * (self._degree + 2) // 2
        if isinstance(self._cell, Quadrilateral):
            return (self._degree + 1) ** 2
        if isinstance(self._cell, Tetrahedron):
            return (self._degree + 1) * (self._degree + 2) * (self._degree + 3) // 6
        if isinstance(self._cell, Hexahedron):
            return (self._degree + 1) ** 3
        raise RuntimeError("Unsupported cell type")

    @property
    def map_type(self) -> AbstractReferenceMap:
        """Get the push forward and pull back map."""
        return IdentityReferenceMap()

    def __hash__(self):
        """Hash."""
        return hash(("uflx_test.LagrangeElement", self._cell, self._degree))

    def tabulate(self, points: np.ndarray, derivative: tuple[int, ...]) -> np.ndarray:
        """Create table of basis function values."""
        if isinstance(self._cell, Point):
            return np.array([[1.0] for () in points])

        if isinstance(self._cell, Interval):
            if self._degree == 0:
                if derivative == (0,):
                    return np.array([[1] for (x,) in points])
                return np.array([[0] for (x,) in points])
            if self._degree == 1:
                if derivative == (0,):
                    return np.array([[1 - x, x] for (x,) in points])
                if derivative == (1,):
                    return np.array([[-1, 1] for (x,) in points])
                return np.array([[0, 0] for (x,) in points])
            if self._degree == 2:
                if derivative == (0,):
                    return np.array(
                        [
                            [(2 * x - 1) * (x - 1), x * (2 * x - 1), 4 * x * (1 - x)]
                            for (x,) in points
                        ]
                    )
                if derivative == (1,):
                    return np.array([[4 * x - 3, 4 * x - 1, 4 - 8 * x] for (x,) in points])

        if isinstance(self._cell, Triangle):
            if self._degree == 0:
                if derivative == (0, 0):
                    return np.array([[1] for (x, y) in points])
                return np.array([[0] for (x, y) in points])
            if self._degree == 1:
                if derivative == (0, 0):
                    return np.array([[1 - x - y, x, y] for (x, y) in points])
                if derivative == (1, 0):
                    return np.array([[-1, 1, 0] for (x, y) in points])
                if derivative == (0, 1):
                    return np.array([[-1, 0, 1] for (x, y) in points])
            if self._degree == 2:
                if derivative == (0, 0):
                    return np.array(
                        [
                            [
                                (1 - x - y) * (1 - 2 * x - 2 * y),
                                x * (2 * x - 1),
                                y * (2 * y - 1),
                                4 * x * y,
                                4 * y * (1 - x - y),
                                4 * x * (1 - x - y),
                            ]
                            for (x, y) in points
                        ]
                    )
                if derivative == (1, 0):
                    return np.array(
                        [
                            [-3 + 4 * x + 4 * y, 4 * x - 1, 0, 4 * y, -4 * y, 4 - 8 * x - 4 * y]
                            for (x, y) in points
                        ]
                    )
                if derivative == (0, 1):
                    return np.array(
                        [
                            [-3 + 4 * x + 4 * y, 0, 4 * y - 1, 4 * x, 4 - 4 * x - 8 * y, -4 * x]
                            for (x, y) in points
                        ]
                    )

        if isinstance(self._cell, Quadrilateral):
            if self._degree == 0:
                if derivative == (0, 0):
                    return np.array([[1] for (x, y) in points])
                return np.array([[0] for (x, y) in points])
            if self._degree == 1:
                if derivative == (0, 0):
                    return np.array(
                        [[(1 - x) * (1 - y), x * (1 - y), (1 - x) * y, x * y] for (x, y) in points]
                    )
                if derivative == (1, 0):
                    return np.array([[y - 1, 1 - y, -y, y] for (x, y) in points])
                if derivative == (0, 1):
                    return np.array([[x - 1, -x, 1 - x, x] for (x, y) in points])

        raise NotImplementedError


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

    def __hash__(self):
        """Hash."""
        return hash("uflx_test.Point")


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

    def __hash__(self):
        """Hash."""
        return hash("uflx_test.Interval")


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

    def __hash__(self):
        """Hash."""
        return hash("uflx_test.Triangle")


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

    def __hash__(self):
        """Hash."""
        return hash("uflx_test.Quadrilateral")


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

    def __hash__(self):
        """Hash."""
        return hash("uflx_test.Tetrahedron")


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

    def __hash__(self):
        """Hash."""
        return hash("uflx_test.Hexahedron")


point = Point()
interval = Interval()
triangle = Triangle()
tetrahedron = Tetrahedron()
quadrilateral = Quadrilateral()
hexahedron = Hexahedron()
