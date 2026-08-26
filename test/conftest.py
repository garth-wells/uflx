"""Implmentations of domains."""

import pytest

from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractReferenceMappedFiniteElement
from uflx.maps import AbstractReferenceMap, IdentityReferenceMap


class Entity(AbstractEntity):
    """An entity."""

    def __init__(self, sub_entities: list[list[tuple[AbstractEntity, tuple[int, ...]]]]):
        """Initalise."""
        self._sub_entities = sub_entities

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

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""
        return len(self._sub_entities) - 1

    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension."""
        return [i[0] for i in self._sub_entities[dim]]

    def sub_entity_vertices(self, dim: int) -> list[list[int]]:
        """Get lists of the vertices of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            A list of lists of vertices of sub-entities of the given dimension.
        """
        return [list(i[1]) for i in self._sub_entities[dim]]


class Point(Entity):
    """A point."""

    def __init__(self):
        """Initalise."""
        super().__init__([[(self, (0,))]])


class Interval(Entity):
    """An interval."""

    def __init__(self):
        """Initalise."""
        super().__init__([[(Point(), (i,)) for i in range(2)], [(self, (0,))]])


class Triangle(Entity):
    """A triangle cell."""

    def __init__(self):
        """Initalise."""
        super().__init__(
            [
                [(Point(), (i,)) for i in range(3)],
                [(Interval(), vs) for vs in [(0, 1), (0, 2), (1, 2)]],
                [(self, (0,))],
            ]
        )


class Quadrilateral(Entity):
    """A quadrilateral."""

    def __init__(self):
        """Initalise."""
        super().__init__(
            [
                [(Point(), (i,)) for i in range(4)],
                [(Interval(), vs) for vs in [(0, 1), (0, 2), (1, 3), (2, 3)]],
                [(self, (0,))],
            ]
        )


class Tetrahedron(Entity):
    """A tetrahedron."""

    def __init__(self):
        """Initalise."""
        super().__init__(
            [
                [(Point(), (i,)) for i in range(4)],
                [(Interval(), vs) for vs in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]],
                [(Triangle(), vs) for vs in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]],
                [(self, (0,))],
            ]
        )


class Hexahedron(Entity):
    """A hexahedron."""

    def __init__(self):
        """Initalise."""
        super().__init__(
            [
                [(Point(), (i,)) for i in range(8)],
                [
                    (Interval(), vs)
                    for vs in [
                        (0, 1),
                        (0, 2),
                        (0, 4),
                        (1, 3),
                        (1, 5),
                        (2, 3),
                        (2, 6),
                        (3, 7),
                        (4, 5),
                        (4, 6),
                        (5, 7),
                        (6, 7),
                    ]
                ],
                [
                    (Quadrilateral(), vs)
                    for vs in [
                        (0, 1, 2, 3),
                        (0, 1, 4, 5),
                        (0, 2, 4, 6),
                        (1, 3, 5, 7),
                        (2, 3, 6, 7),
                        (4, 5, 6, 7),
                    ]
                ],
                [(self, (0,))],
            ]
        )


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
    """Create a Lagrange element."""

    def create(cell_name: str, degree: int, block_shape: tuple[int, ...] | None = None):
        match cell_name:
            case "point":
                cell: AbstractEntity = Point()
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
    """Generate tests."""
    if "entity" in metafunc.fixturenames:
        metafunc.parametrize(
            "entity",
            [
                pytest.param(Point(), id="point"),
                pytest.param(Interval(), id="interval"),
                pytest.param(Triangle(), id="triangle"),
                pytest.param(Quadrilateral(), id="quadrilateral"),
                pytest.param(Tetrahedron(), id="tetrahedron"),
                pytest.param(Hexahedron(), id="hexahedron"),
            ],
        )
    if "simplex" in metafunc.fixturenames:
        metafunc.parametrize(
            "simplex",
            [
                pytest.param(Point(), id="point"),
                pytest.param(Interval(), id="interval"),
                pytest.param(Triangle(), id="triangle"),
                pytest.param(Tetrahedron(), id="tetrahedron"),
            ],
        )
