"""Implmentations of domains."""

import pytest
import numpy as np

from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractReferenceMappedFiniteElement
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


class Quadrilateral(AbstractEntity):
    """A quadrilateral."""

    @property
    def name(self) -> str:
        """Name of this entity."""
        return "quadrilateral"

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
        return hash("uflx.test.Quadrilateral")


class Tetrahedron(AbstractEntity):
    """A tetrahedron."""

    @property
    def name(self) -> str:
        """Name of this entity."""
        return "tetrahedron"

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
        return hash("uflx.test.Tetrahedron")


class Hexahedron(AbstractEntity):
    """A hexahedron."""

    @property
    def name(self) -> str:
        """Name of this entity."""
        return "hexahedron"

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
                return [Point() for _ in range(8)]
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
        return hash("uflx.test.Hexahedron")


class LagrangeElement(AbstractReferenceMappedFiniteElement):
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
        if isinstance(self._cell, Quadrilateral):
            return (self._degree + 1) ** 2
        if isinstance(self._cell, Tetrahedron):
            return (self._degree + 1) * (self._degree + 2) * (self._degree + 3) // 6
        if isinstance(self._cell, Hexahedron):
            return (self._degree + 1) ** 3
        raise RuntimeError("Unsupported cell type")

    @property
    def reference_map(self) -> AbstractReferenceMap:
        """Get the push forward and pull back map."""
        return IdentityReferenceMap()

    def __hash__(self):
        """Hash."""
        return hash(("uflx_test.LagrangeElement", self._cell, self._degree))


@pytest.fixture
def lagrange_element():
    def create(cell_name: str, degree: int, block_shape: tuple[int, ...] | None = None):
        match cell_name:
            case "point":
                cell = Point()
            case "interval":
                cell = Interval()
            case "triangle":
                cell = Triangle()
            case "quadrilateral":
                cell = Quadrilateral()
            case "tetrahedron":
                cell = Tetrahedron()
            case "hexahedron":
                cell = Hexahedron()
            case _:
                raise ValueError(f"Invalid cell: {cell_name}")
        return LagrangeElement(cell, degree, block_shape)
    return create


def pytest_generate_tests(metafunc):
    if "entity" in metafunc.fixturenames:
        metafunc.parametrize("entity", [
            pytest.param(Point(), id="point"),
            pytest.param(Interval(), id="interval"),
            pytest.param(Triangle(), id="triangle"),
            pytest.param(Quadrilateral(), id="quadrilateral"),
            pytest.param(Tetrahedron(), id="tetrahedron"),
            pytest.param(Hexahedron(), id="hexahedron"),
        ])
    if "simplex" in metafunc.fixturenames:
        metafunc.parametrize("simplex", [
            pytest.param(Point(), id="point"),
            pytest.param(Interval(), id="interval"),
            pytest.param(Triangle(), id="triangle"),
            pytest.param(Tetrahedron(), id="tetrahedron"),
        ])
