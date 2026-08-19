"""Implmentations of finite elements."""

import numpy as np

from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractReferenceMappedFiniteElement
from uflx.maps import AbstractReferenceMap, IdentityReferenceMap
from uflx.test_utils.domains import (
    Hexahedron,
    Interval,
    Point,
    Quadrilateral,
    Tetrahedron,
    Triangle,
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
