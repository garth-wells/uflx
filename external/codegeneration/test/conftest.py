"""Implmentations of domains."""

import numpy as np
import pytest
from uflx.entities import AbstractEntity
from uflx_codegeneration.finite_element import FiniteElement
from uflx.maps import AbstractReferenceMap, IdentityReferenceMap


class Point(AbstractEntity):
    """A point."""

    @property
    def name(self) -> str:
        """Name of this entity."""
        return "point"

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
        return hash("uflx.test.Point")


class Interval(AbstractEntity):
    """An interval."""

    @property
    def name(self) -> str:
        """Name of this entity."""
        return "interval"

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
        return hash("uflx.test.Interval")


class Triangle(AbstractEntity):
    """A triangle cell."""

    @property
    def name(self) -> str:
        """Name of this entity."""
        return "triangle"

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
        return hash("uflx.test.Triangle")


class LagrangeElement(FiniteElement):
    """A Lagrange element."""

    def __init__(
        self, cell: AbstractEntity, degree: int, block_shape: tuple[int, ...] | None = None
    ):
        """Create a Lagrange element."""
        self._cell = cell
        self._degree = degree
        self._block_shape = block_shape

    def __repr__(self):
        """Representation."""
        return f"uflx.test.LagrangeElement({self._cell!r}, {self._degree}, self._block_shape)"

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

    @property
    def lagrange_superdegree(self) -> int | None:
        """Degree of the minimum degree Lagrange space that spans this element."""
        return self._degree

    @property
    def dim(self) -> int:
        """The dimension of the finite element, ie the number of basis functions."""
        if isinstance(self._cell, Point):
            return 1
        if isinstance(self._cell, Interval):
            return self._degree + 1
        if isinstance(self._cell, Triangle):
            return (self._degree + 1) * (self._degree + 2) // 2
        raise RuntimeError("Unsupported cell type")

    @property
    def reference_map(self) -> AbstractReferenceMap:
        """Get the push forward and pull back map."""
        return IdentityReferenceMap()

    def __hash__(self):
        """Hash."""
        return hash(("uflx_test.LagrangeElement", self._cell, self._degree))

    def tabulate(self, points: np.ndarray, derivative: tuple[int, ...]) -> np.ndarray:
        """Create table of basis function values."""
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

        raise NotImplementedError


@pytest.fixture
def lagrange_element():
    """Create a Lagrange element."""

    def create(cell_name: str, degree: int, block_shape: tuple[int, ...] | None = None):
        match cell_name:
            case "point":
                cell: AbstractEntity = Point()
            case "interval":
                cell = Interval()
            case "triangle":
                cell = Triangle()
            case _:
                raise ValueError(f"Invalid cell: {cell_name}")
        return LagrangeElement(cell, degree, block_shape)

    return create
