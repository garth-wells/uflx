"""Implmentations of domains."""

import os

import numpy as np
import numpy.typing as npt
import pytest
from uflx.entities import AbstractEntity
from uflx.maps import AbstractReferenceMap, IdentityReferenceMap

from uflx_codegeneration.finite_element import AbstractFiniteElement, BlockedElement
from uflx_codegeneration.utils import number_of_derivatives


class Entity(AbstractEntity):
    """An entity."""

    def __repr__(self):
        """Representation."""
        return self.__class__.__name__

    def __str__(self):
        """String."""
        return f"{self!r}"

    def __hash__(self):
        """Hash."""
        return hash(("uflx.test", self.__class__.__name__))

    @property
    def name(self) -> str:
        """Name of this entity."""
        return self.__class__.__name__.lower()

    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""
        return isinstance(other, Entity) and f"{self!r}" == f"{other!r}"


class Point(Entity):
    """A point."""

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


class Interval(Entity):
    """An interval."""

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


class Triangle(Entity):
    """A triangle cell."""

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


class LagrangeElement(AbstractFiniteElement):
    """A Lagrange element."""

    def __init__(self, cell: AbstractEntity, degree: int):
        """Create a Lagrange element."""
        self._cell = cell
        self._degree = degree

    def __repr__(self):
        """Representation."""
        return f"uflx.test.LagrangeElement({self._cell!r}, {self._degree})"

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

    @property
    def lagrange_superdegree(self) -> int:
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

    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        """Create table of basis function values and derivatives."""
        points = np.asarray(points)
        table = np.zeros(
            [number_of_derivatives(derivatives, self.cell), points.shape[0], self.dim, 1]
        )

        if isinstance(self._cell, Triangle):
            match self._degree:
                case 0:
                    table[0, :, :, :] = 1.0
                case 1:
                    table[0, :, :, 0] = [[1 - x - y, x, y] for (x, y) in points]
                    if derivatives >= 1:
                        table[1, :, :, 0] = [[-1, 1, 0] for (x, y) in points]
                        table[2, :, :, 0] = [[-1, 0, 1] for (x, y) in points]
                case 2:
                    table[0, :, :, 0] = [
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
                    if derivatives >= 1:
                        table[1, :, :, 0] = [
                            [-3 + 4 * x + 4 * y, 4 * x - 1, 0, 4 * y, -4 * y, 4 - 8 * x - 4 * y]
                            for (x, y) in points
                        ]
                        table[2, :, :, 0] = [
                            [-3 + 4 * x + 4 * y, 0, 4 * y - 1, 4 * x, 4 - 4 * x - 8 * y, -4 * x]
                            for (x, y) in points
                        ]
                    if derivatives >= 2:
                        table[3, :, :, 0] = [[4, 4, 0, 0, 0, -8 * x] for (x, y) in points]
                        table[4, :, :, 0] = [[4, 0, 0, 4, -4, -4] for (x, y) in points]
                        table[5, :, :, 0] = [[4, 0, 4, 0, -8, 0] for (x, y) in points]
                case _:
                    raise NotImplementedError()
            return table

        raise NotImplementedError()


@pytest.fixture
def lagrange_element():
    """Create a Lagrange element."""

    def create(cell_name: str, degree: int, block_shape: tuple[int, ...] | None = None):
        if block_shape is not None:
            return BlockedElement(create(cell_name, degree), block_shape)
        match cell_name:
            case "point":
                cell: AbstractEntity = Point()
            case "interval":
                cell = Interval()
            case "triangle":
                cell = Triangle()
            case _:
                raise ValueError(f"Invalid cell: {cell_name}")
        return LagrangeElement(cell, degree)

    return create


CODE_DIRECTORY = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), ".code")


@pytest.fixture
def code_dir():
    """Directory where code generated by tests is saved."""
    return CODE_DIRECTORY


def pytest_generate_tests(metafunc):
    """Generate tests."""
    if not os.path.isdir(CODE_DIRECTORY):
        os.mkdir(CODE_DIRECTORY)
