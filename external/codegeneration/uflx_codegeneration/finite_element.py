"""Finite element."""

from abc import abstractmethod
from collections.abc import Sequence

import numpy.typing as npt
from uflx.finite_elements import AbstractReferenceMappedFiniteElement


def product(ls: Sequence[int]) -> int:
    """Return the product of numbers in a list."""
    result = 1
    for i in ls:
        result *= i
    return result


class FiniteElement(AbstractReferenceMappedFiniteElement):
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
    def reference_value_size(self) -> int:
        """Return the value size of the value space on the reference cell."""
        return product(self.reference_value_shape)

    def physical_value_size(self, geometric_dimension: int) -> int:
        """Return the value size of the value space on a physical cell."""
        return product(self.physical_value_shape(geometric_dimension))

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
