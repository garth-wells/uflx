"""Utilities for testing UFL."""

import typing

from uflx.finite_element import AbstractFiniteElement
from uflx.cell import AbstractCell


class LagrangeElement(AbstractFiniteElement):
    """A Lagrange element."""

    def __init__(self, cell: AbstractCell, degree: int):
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
    def cell(self) -> AbstractCell:
        """Return the cell that this element is defined on."""
        return self._cell

    @property
    def reference_value_shape(self) -> typing.Tuple[int, ...]:
        """Return the shape of the value space on the reference cell."""
        return ()

    @property
    def embedded_lagrange_superdegree(self) -> int | None:
        """Degree of the minimum degree Lagrange space that spans this element."""
        return self._degree


class Point(AbstractCell):
    """A point."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Point)

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return 0

    def sub_entities(self, dim: int) -> typing.List[AbstractCell]:
        """Get a list of sub-entities of a given dimension."""
        match dim:
            case 0:
                return [self]
            case _:
                raise ValueError(f"Invalid dimension: {dim}")


class Interval(AbstractCell):
    """An interval."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Interval)

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return 1

    def sub_entities(self, dim: int) -> typing.List[AbstractCell]:
        """Get a list of sub-entities of a given dimension."""
        match dim:
            case 0:
                return [Point() for _ in range(2)]
            case 1:
                return [self]
            case _:
                raise ValueError(f"Invalid dimension: {dim}")


class Triangle(AbstractCell):
    """A triangle."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Triangle)

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return 2

    def sub_entities(self, dim: int) -> typing.List[AbstractCell]:
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


class Tetrahedron(AbstractCell):
    """A tetrahedron."""

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Tetrahedron)

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return 3

    def sub_entities(self, dim: int) -> typing.List[AbstractCell]:
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


point = Point()
interval = Interval()
triangle = Triangle()
tetrahedron = Tetrahedron()
