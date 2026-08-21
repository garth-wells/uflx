"""Finite element."""

from abc import abstractmethod
from enum import Enum

import numpy as np
import numpy.typing as npt
from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractReferenceMappedFiniteElement
from uflx.maps import (
    AbstractReferenceMap,
    BlockedReferenceMap,
    MixedReferenceMap,
    SymmetricReferenceMap,
)


class AbstractFiniteElement(AbstractReferenceMappedFiniteElement):
    """A base wrapper for a uflx-codegeneration compatible finite element.

    This class includes methods and properties required for code generation.
    By implementing these functions in a finite element definition, you can
    enable that finite element library to be used to generate finite element
    kernel code.

    Note that this class inherits from UFLx's AbstractReferenceMappedFiniteElement,
    and so all the abstract methods from that class must be implemented too.
    The return type of the property `lagrange_superdegree` that is defined in
    this class differs from the return type of AbstractReferenceMappedFiniteElement:
    this property cannot be None in order for code to be successfully generated.
    """

    @abstractmethod
    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        """Create table of basis function values and derivatives.

        This function returns a four dimensional array whose shape is (number
        of derivatives, number of points, number of basis functions, value size).
        The derivatives are sotred in triangular (2D) or tetrahedral (3D)
        ordering - for example, in 2D the derivatives are in the following order,
        where ``(x,y)`` represents ``x`` derivatives in the x-direction and ``y``
        in the y-direction:
        ``(0,0)``, ``(1,0)``, ``(0,1)``, ``(2,0)``, ``(1,1)``, ``(0,2)``, ``(3,0)``,
        ...
        """

    @property
    @abstractmethod
    def lagrange_superdegree(self) -> int:
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


class MixedElement(AbstractFiniteElement):
    """A mixed element that combines two or more elements.

    This can be used when multiple different elements appear in a form.
    Users should not directly call this class's initilizer, but should
    use the :func:`mixed_element` function instead.
    """

    _sub_elements: list[AbstractFiniteElement]

    def __init__(self, sub_elements: list[AbstractFiniteElement]):
        """Initialise the element."""
        assert len(sub_elements) > 0
        self._sub_elements = sub_elements
        self._reference_map = MixedReferenceMap(
            [e.reference_map for e in sub_elements], [e.reference_value_shape for e in sub_elements]
        )
        for e in sub_elements[1:]:
            if e.cell != sub_elements[0].cell:
                raise ValueError(
                    "Cannot create mixed element where sub-elements are defined on different cells."
                )

    def __hash__(self):
        """Hash."""
        return hash(("uflx_codegeneration", f"{self!r}"))

    def __repr__(self) -> str:
        """Representation."""
        return "MixedElement(" + ", ".join(f"{i!r}" for i in self._sub_elements) + ")"

    def __str__(self):
        """String."""
        return f"{self!r}"

    def __eq__(self, other) -> bool:
        """Check if this element is equal to another element."""
        return (
            isinstance(other, MixedElement)
            and len(self._sub_elements) == len(other._sub_elements)
            and all(i == j for i, j in zip(self._sub_elements, other._sub_elements))
        )

    @property
    def cell(self) -> AbstractEntity:
        """Return the cell that this element is defined on."""
        return self._sub_elements[0].cell

    @property
    def dim(self) -> int:
        """The dimension of the finite element, ie the number of basis functions."""
        return sum(e.dim for e in self._sub_elements)

    @property
    def lagrange_superdegree(self) -> int:
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
        return max([e.lagrange_superdegree for e in self._sub_elements], default=0)

    @property
    def reference_map(self) -> AbstractReferenceMap:
        """Get the push forward and pull back map."""
        return self._reference_map

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        """Return the shape of the value space on the reference cell."""
        raise NotImplementedError()

    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        """Create table of basis function values and derivatives.

        This function returns a four dimensional array whose shape is (number
        of derivatives, number of points, number of basis functions, value size).
        The derivatives are sotred in triangular (2D) or tetrahedral (3D)
        ordering - for example, in 2D the derivatives are in the following order,
        where ``(x,y)`` represents ``x`` derivatives in the x-direction and ``y``
        in the y-direction:
        ``(0,0)``, ``(1,0)``, ``(0,1)``, ``(2,0)``, ``(1,1)``, ``(0,2)``, ``(3,0)``,
        ...
        """
        return np.stack([e.tabulate(derivatives, points) for e in self._sub_elements], axis=2)


class BlockedOrdering(Enum):
    """Ordering of components in a blocked element."""

    xyzxyz = 1
    xxyyzz = 2


class BlockedElement(AbstractFiniteElement):
    """Element with a block size that contains multiple copies of a sub element.

    This can be used to (for example) create vector and tensor Lagrange
    elements. Users should not directly call this classes initilizer,
    but should use the `blocked_element` function instead.

    """

    _block_shape: tuple[int, ...]
    _sub_element: AbstractFiniteElement
    _block_size: int
    _has_symmetry: bool

    def __init__(
        self,
        sub_element: AbstractFiniteElement,
        shape: tuple[int, ...],
        symmetry: bool | None = None,
        ordering: BlockedOrdering = BlockedOrdering.xyzxyz,
    ):
        """Initialise the element."""
        if sub_element.reference_value_shape != ():
            raise ValueError(
                "Blocked elements of non-scalar elements are not supported. "
                "Try using MixedElement instead."
            )
        if symmetry is not None:
            if len(shape) != 2:
                raise ValueError("symmetry argument can only be passed to elements of rank 2.")
            if shape[0] != shape[1]:
                raise ValueError("symmetry argument can only be passed to square shaped elements.")

        if symmetry:
            block_size = shape[0] * (shape[0] + 1) // 2
            self._has_symmetry = True
        else:
            block_size = 1
            for i in shape:
                block_size *= i
            self._has_symmetry = False
        assert block_size > 0

        self._sub_element = sub_element
        self._block_size = block_size
        self._block_shape = shape

        if symmetry:
            n = 0
            symmetry_map: dict[tuple[int, ...], int] = {}
            for i in range(shape[0]):
                for j in range(i + 1):
                    symmetry_map[(i, j)] = n
                    symmetry_map[(j, i)] = n
                    n += 1

            self._reference_map = SymmetricReferenceMap(
                sub_element.reference_map,
                shape,
                symmetry_map,
            )
        self._symmetry = symmetry
        self._ordering = ordering

    @property
    def sub_element(self):
        """Get the scalar sub-element."""
        return self._sub_element

    def __hash__(self):
        """Hash."""
        return hash(("uflx_codegeneration", f"{self!r}"))

    def __repr__(self) -> str:
        """Representation."""
        repr = f"BlockedElement({self._sub_element!r}, {self._block_shape}"
        if self._symmetry is not None:
            repr += f", symmetry={self._symmetry}"
        repr += ")"
        return repr

    def __str__(self):
        """String."""
        return f"{self!r}"

    def __eq__(self, other) -> bool:
        """Check if this element is equal to another element."""
        return (
            isinstance(other, BlockedElement)
            and self._block_size == other._block_size
            and self._block_shape == other._block_shape
            and self._sub_element == other._sub_element
        )

    @property
    def cell(self) -> AbstractEntity:
        """Return the cell that this element is defined on."""
        return self._sub_element.cell

    @property
    def dim(self) -> int:
        """The dimension of the finite element, ie the number of basis functions."""
        return self._sub_element.dim * self._block_size

    @property
    def lagrange_superdegree(self) -> int:
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
        return self._sub_element.lagrange_superdegree

    @property
    def reference_map(self) -> AbstractReferenceMap:
        """Get the push forward and pull back map."""
        return BlockedReferenceMap(self._sub_element.reference_map, self._block_shape)

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        """Return the value size of the value space on the reference cell."""
        return self._block_shape

    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        """Create table of basis function values and derivatives.

        This function returns a four dimensional array whose shape is (number
        of derivatives, number of points, number of basis functions, value size).
        The derivatives are sotred in triangular (2D) or tetrahedral (3D)
        ordering - for example, in 2D the derivatives are in the following order,
        where ``(x,y)`` represents ``x`` derivatives in the x-direction and ``y``
        in the y-direction:
        ``(0,0)``, ``(1,0)``, ``(0,1)``, ``(2,0)``, ``(1,1)``, ``(0,2)``, ``(3,0)``,
        ...
        """
        scalar_table = self._sub_element.tabulate(derivatives, points)
        table = np.zeros(
            (
                scalar_table.shape[0],
                scalar_table.shape[1],
                scalar_table.shape[2] * self._block_size,
                scalar_table.shape[3] * self._block_size,
            )
        )
        match self._ordering:
            case BlockedOrdering.xyzxyz:
                d = scalar_table.shape[3]
                for i in range(self._block_size):
                    table[:, :, i :: self._block_size, i * d : (i + 1) * d] = scalar_table
            case BlockedOrdering.xxyyzz:
                d = scalar_table.shape[2]
                for i in range(self._block_size):
                    table[:, :, i * d : (i + 1) * d, i :: self._block_size] = scalar_table
        return table
