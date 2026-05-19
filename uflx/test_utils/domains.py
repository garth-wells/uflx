"""Implmentations of domains"""

from uflx.entities import AbstractEntity


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
