"""Basix UFLx interface."""

from abc import abstractmethod
from collections.abc import Sequence

import hashlib
from warnings import warn

import numpy as np
import numpy.typing as npt
import uflx
from uflx.finite_elements import AbstractReferenceMappedFiniteElement as ARMFE
from uflx.scalars import AbstractInteger
from uflx_codegeneration import FiniteElement

import basix

__all__ = [
    "element",
    "custom_element",
    "mixed_element",
    "quadrature_element",
    "real_element",
    "blocked_element",
    "wrap_element",
]


def convert_map(basix_map: basix.MapType) -> uflx.maps.AbstractReferenceMap:
    """Convert a basix map to a UFLx map."""
    match basix_map:
        case basix.MapType.identity:
            return uflx.maps.IdentityReferenceMap()
        case _:
            raise NotImplementedError()


class BasixCell(uflx.entities.AbstractEntity):
    """A cell defined by Basix."""

    def __init__(self, basix_cell: basix.CellType):
        self._basix_cell = basix_cell

    def __eq__(self, other) -> bool:
        """Check if this entity is equal to another entity."""
        return isinstance(other, BasixCell) and self._basix_cell == other.basix_cell

    @property
    def topological_dimension(self) -> int:
        """Topological dimension of the entity."""
        return len(basix.cell.topology(self._basix_cell)) - 1

    def sub_entities(self, dim: int) -> list[uflx.entities.AbstractEntity]:
        """Get a list of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            A list of sub-entities of the given dimension.
        """
        return basix.cell.subentity_types(self._basix_cell)[dim]

    def __hash__(self):
        return hash(("basix.uflx", f"{self!r}"))

    def __repr__(self):
        return f"BasixCell({self._basix_cell.name})"

    def __str__(self):
        return f"{self!r}"


def hash_data(data: Sequence[np.floating] | np.floating):
    """Return a hash of an array of floating point numbers"""

    def hash_data_inner(data: Sequence[np.floating] | np.floating):
        """Represent an array as a string, ready for hashing."""
        try:
            return ",".join(hash_data_inner(i) for i in data)  # type: ignore
        except TypeError:
            return f"{data!r}"

    return hashlib.sha1(hash_data_inner(data).encode("utf-8")).hexdigest()


class BasixElement(FiniteElement):
    """A wrapper allowing Basix elements to be used directly with UFLx.

    This class allows elements created with `basix.create_element` to be
    wrapped as UFLx compatible elements. Users should not directly call
    this class's initialiser, but should use the `element` function
    instead.
    """

    _element: basix.finite_element.FiniteElement

    def __init__(self, element: basix.finite_element.FiniteElement):
        """Initialise."""
        self._element = element

    def __hash__(self):
        return hash(("basix.uflx", f"{self!r}"))

    def __repr__(self) -> str:
        if self._element.family == basix.ElementFamily.custom:
            return (
                f"FiniteElement(CUSTOM, {self._element.cell_type.name}, "
                f"{self._element.value_shape}, {self._element.map_type.name}, "
                f"{self._element.discontinuous}, {self._element.embedded_subdegree}, "
                f"{self._element.embedded_superdegree}, {self._element.dtype}, "
                f"{self._element.dof_ordering}, {hash_data(self._element.wcoeffs)}, "
                f"{hash_data(self._element.x)}, {hash_data(self._element.M)}"
            )
        else:
            return (
                f"FiniteElement({self._element.family.name}, {self._element.cell_type.name}, "
                f"{self._element.degree}, "
                f"{self._element.lagrange_variant.name}, {self._element.dpc_variant.name}, "
                f"{self._element.discontinuous}, "
                f"{self._element.dtype}, {self._element.dof_ordering})"
            )

    def __str__(self):
        return f"{self!r}"

    def __eq__(self, other) -> bool:
        return isinstance(other, BasixElement) and self._element == other._element

    @property
    def cell(self) -> uflx.entities.AbstractEntity:
        return BasixCell(self._element.cell_type)

    @property
    def dim(self) -> int | AbstractInteger:
        return self._element.dim

    @property
    def lagrange_superdegree(self) -> int:
        return self._element.embedded_superdegree

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        raise NotImplementedError()

    @property
    def reference_map(self) -> uflx.maps.AbstractReferenceMap:
        return convert_map(self._element.map_type)

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        return tuple(self._element.value_shape)

    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        return self._element.tabulate(derivatives, points)


class MixedElement(FiniteElement):
    """A mixed element that combines two or more elements.

    This can be used when multiple different elements appear in a form.
    Users should not directly call this class's initilizer, but should
    use the :func:`mixed_element` function instead.
    """

    _sub_elements: list[FiniteElement]

    def __init__(self, sub_elements: list[FiniteElement]):
        """Initialise the element."""
        assert len(sub_elements) > 0
        self._sub_elements = sub_elements
        reference_map = (
            uflx.identity_pullback
            if all(isinstance(e.pullback, _IdentityPullback) for e in sub_elements)
            else _MixedPullback(self)
        )

    def __hash__(self):
        return hash(("basix.uflx", f"{self!r}"))

    def __repr__(self) -> str:
        return "MixedElement(" + ", ".join(f"{i!r}" for i in self._sub_elements) + ")"

    def __str__(self):
        return f"{self!r}"


class BlockedElement(FiniteElement):
    """Element with a block size that contains multiple copies of a sub element.

    This can be used to (for example) create vector and tensor Lagrange
    elements. Users should not directly call this classes initilizer,
    but should use the `blocked_element` function instead.

    """

    _block_shape: tuple[int, ...]
    _sub_element: FiniteElement
    _block_size: int
    _has_symmetry: bool

    def __init__(
        self,
        sub_element: FiniteElement,
        shape: tuple[int, ...],
        symmetry: bool | None = None,
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

            self._reference_map = uflx.maps.SymmetricReferenceMap(
                sub_element.reference_map, shape, symmetry_map,
            )
        self._symmetry = symmetry

    def __hash__(self):
        return hash(("basix.uflx", f"{self!r}"))

    def __repr__(self) -> str:
        repr = f"BlockedElement({self._sub_element!r}, {self._block_shape}"
        if self._symmetry is not None:
            repr += f", symmetry={self._symmetry}"
        repr += ")"
        return repr

    def __str__(self):
        return f"{self!r}"

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, BlockedElement)
            and self._block_size == other._block_size
            and self._block_shape == other._block_shape
            and self._sub_element == other._sub_element
        )

    @property
    def cell(self) -> uflx.entities.AbstractEntity:
        return BasixCell(self._sub_element.cell_type)

    @property
    def dim(self) -> int | AbstractInteger:
        return self._sub_element.dim * self._block_size

    @property
    def lagrange_superdegree(self) -> int:
        return self._sub_element.embedded_superdegree

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        return self._block_shape

    @property
    def reference_map(self) -> uflx.maps.AbstractReferenceMap:
        return convert_map(self._sub_element.map_type)

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        return self._block_shape

    def tabulate(self, derivatives: int, points: npt.ArrayLike) -> npt.NDArray:
        scalar_table = self._sub_element.tabulate(derivatives, points)
        table = np.zeros((
            scalar_table.shape[0],
            scalar_table.shape[1],
            scalar_table.shape[2] * self._block_size,
            scalar_table.shape[3] * self._block_size,
        ))
        d = scalar_table.shape[3]
        for i in range(self._block_size):
            table[:, :, i::self._block_size, d*i:d*(i+1)] = scalar_table
        return table



class QuadratureElement(FiniteElement):
    """A quadrature element."""

    def __init__(
        self,
        cell: basix.CellType,
        points: npt.NDArray[np.floating],
        weights: npt.NDArray[np.floating],
        reference_map: uflx.maps.AbstractReferenceMap,
        degree: int | None = None,
        dtype: npt.DTypeLike = np.float64,
    ):
        """Initialise the element."""
        self._points = points.astype(dtype)
        self._weights = weights.astype(dtype)
        self._cell_type = cell
        self._entity_counts = [len(i) for i in basix.topology(cell)]
        self._reference_map = reference_map

        if degree is None:
            degree = len(points)

    def __hash__(self):
        return hash(("basix.uflx", f"{self!r}"))

    def __repr__(self) -> str:
        return (
            f"QuadratureElement({cell.name}, {hash_data(points)}, "
            f"{hash_data(weights)}, {reference_map!r})"
        )

    def __str__(self):
        return f"{self!r}"



class RealElement(FiniteElement):
    """A real element."""

    def __init__(self, cell: basix.CellType, value_shape: tuple[int, ...]):
        """Initialise the element."""
        self._cell_type = cell
        tdim = len(basix.topology(cell)) - 1
        self._value_shape = value_shape

        self._entity_counts = []
        if tdim >= 1:
            self._entity_counts.append(self.cell.num_vertices)
        if tdim >= 2:
            self._entity_counts.append(self.cell.num_edges)
        if tdim >= 3:
            self._entity_counts.append(self.cell.num_facets)
        self._entity_counts.append(1)

    def __hash__(self):
        return hash(("basix.uflx", f"{self!r}"))

    def __repr__(self) -> str:
        return f"RealElement({self._cell_type.name}, {self._value_shape})"

    def __str__(self):
        return f"{self!r}"


def element(
    family: basix.ElementFamily | str,
    cell: basix.CellType | str,
    degree: int,
    lagrange_variant: basix.LagrangeVariant = basix.LagrangeVariant.unset,
    dpc_variant: basix.DPCVariant = basix.DPCVariant.unset,
    discontinuous: bool = False,
    shape: tuple[int, ...] | None = None,
    symmetry: bool | None = None,
    dof_ordering: list[int] | None = None,
    dtype: npt.DTypeLike | None = None,
) -> FiniteElement:
    """Create a UFLx compatible element using Basix.

    Args:
        family: Element family/type.
        cell: Element cell type.
        degree: Degree of the finite element.
        lagrange_variant: Variant of Lagrange to be used.
        dpc_variant: Variant of DPC to be used.
        discontinuous: If ``True``, the discontinuous version of the
            element is created.
        shape: Value shape of the element. For scalar-valued families,
            this can be used to create vector and tensor elements.
        symmetry: Set to ``True`` if the tensor is symmetric. Valid for
            rank 2 elements only.
        dof_ordering: Ordering of dofs for ``ElementDofLayout``.
        dtype: Floating point data type.

    Returns:
        A finite element.

    """
    # Conversion of string arguments to types
    if isinstance(cell, str):
        cell = basix.CellType[cell]
    if isinstance(family, str):
        if family.startswith("Discontinuous "):
            family = family[14:]
            discontinuous = True
        if family in ["DP", "DG", "DQ"]:
            family = "P"
            discontinuous = True
        if family == "DPC":
            discontinuous = True

        family = basix.finite_element.string_to_family(family, cell.name)

    # Default variant choices
    EF = basix.ElementFamily
    if lagrange_variant == basix.LagrangeVariant.unset:
        if family == EF.P:
            lagrange_variant = basix.LagrangeVariant.gll_warped
        elif family in [EF.RT, EF.N1E]:
            lagrange_variant = basix.LagrangeVariant.legendre
        elif family in [EF.serendipity, EF.BDM, EF.N2E]:
            lagrange_variant = basix.LagrangeVariant.legendre

    if dpc_variant == basix.DPCVariant.unset:
        if family in [EF.serendipity, EF.BDM, EF.N2E]:
            dpc_variant = basix.DPCVariant.legendre
        elif family == EF.DPC:
            dpc_variant = basix.DPCVariant.diagonal_gll

    e = basix.create_element(
        family,
        cell,
        degree,
        lagrange_variant,
        dpc_variant,
        discontinuous,
        dof_ordering=dof_ordering,
        dtype=dtype,
    )
    ufl_e = BasixElement(e)

    if shape is None or shape == tuple(e.value_shape):
        if symmetry is not None:
            raise ValueError("Cannot pass a symmetry argument to this element.")
        return ufl_e
    else:
        return blocked_element(ufl_e, shape=shape, symmetry=symmetry)


def custom_element(
    cell_type: basix.CellType,
    reference_value_shape: Sequence[int],
    wcoeffs: npt.ArrayLike,
    x: Sequence[Sequence[npt.ArrayLike]],
    M: Sequence[Sequence[npt.ArrayLike]],
    interpolation_nderivs: int,
    map_type: basix.MapType,
    sobolev_space: basix.SobolevSpace,
    discontinuous: bool,
    embedded_subdegree: int,
    embedded_superdegree: int,
    polyset_type: basix.PolysetType = basix.PolysetType.standard,
    dtype: npt.DTypeLike | None = None,
) -> FiniteElement:
    """Create a UFLx compatible custom Basix element.

    Args:
        cell_type: The cell type
        reference_value_shape: The reference value shape of the element
        wcoeffs: Matrices for the kth value index containing the
            expansion coefficients defining a polynomial basis spanning
            the polynomial space for this element. Shape is
            ``(dim(finite element polyset), dim(Legenre polynomials))``.
        x: Interpolation points. Indices are ``(tdim, entity index,
            point index, dim)``.
        M: The interpolation matrices. Indices are ``(tdim, entity
            index, dof, vs, point_index, derivative)``.
        interpolation_nderivs: The number of derivatives that need to be
            used during interpolation.
        map_type: The type of map to be used to map values from the
            reference to a cell.
        sobolev_space: Underlying Sobolev space for the element.
        discontinuous: Indicates whether or not this is the
            discontinuous version of the element.
        embedded_subdegree: The highest degree ``n`` such that a
            Lagrange (or vector Lagrange) element of degree ``n`` is a
            subspace of this element.
        embedded_superdegree: The highest degree of a polynomial in this
            element's polyset.
        polyset_type: Polyset type for the element.
        dtype: Floating point data type.

    Returns:
        A custom finite element.
    """
    e = basix.create_custom_element(
        cell_type,
        tuple(reference_value_shape),
        np.array(wcoeffs),
        [[_np.array(j) for j in i] for i in x],
        [[_np.array(j) for j in i] for i in M],
        interpolation_nderivs,
        map_type,
        sobolev_space,
        discontinuous,
        embedded_subdegree,
        embedded_superdegree,
        polyset_type,
        dtype=dtype,
    )
    return BasixElement(e)


def mixed_element(elements: Sequence[FiniteElement]) -> FiniteElement:
    """Create a UFLx compatible mixed element from a list of elements.

    Args:
        elements: The list of elements

    Returns:
        A mixed finite element.
    """
    return MixedElement(list(elements))


def quadrature_element(
    cell: str | basix.CellType,
    value_shape: Sequence[int] = (),
    scheme: str | None = None,
    degree: int | None = None,
    points: npt.NDArray[np.floating] | None = None,
    weights: npt.NDArray[np.floating] | None = None,
    reference_map: uflx.maps.AbstractReferenceMap | None = None,
    symmetry: bool | None = None,
    dtype: npt.DTypeLike | None = None,
) -> FiniteElement:
    """Create a quadrature element.

    When creating this element, either the quadrature scheme and degree
    must be input or the quadrature points and weights must be.

    Args:
        cell: Cell to create the element on.
        value_shape: Value shape of the element.
        scheme: Quadrature scheme.
        degree: Quadrature degree.
        points: Quadrature points.
        weights: Quadrature weights.
        reference_map: Reference map
        symmetry: Set to ``True`` if the tensor is symmetric. Valid for
            rank 2 elements only.
        dtype: Data type of quadrature points and weights

    Returns:
        A 'quadrature' finite element.
    """
    if isinstance(cell, str):
        cell = basix.CellType[cell]

    if reference_map is None:
        refernece_map = uflx.maps.IdentityReferenceMap()

    if points is None:
        assert weights is None
        assert degree is not None
        if scheme is None:
            points, weights = basix.make_quadrature(cell, degree)  # type: ignore
        else:
            points, weights = basix.make_quadrature(  # type: ignore
                cell, degree, rule=basix.quadrature.string_to_type(scheme)
            )

    assert points is not None
    assert weights is not None

    e = QuadratureElement(cell, points, weights, reference_map, degree, dtype=dtype)
    if tuple(value_shape) == ():
        if symmetry is not None:
            raise ValueError("Cannot pass a symmetry argument to this element.")
        return e
    else:
        return blocked_element(e, shape=tuple(value_shape), symmetry=symmetry)


def real_element(
    cell: basix.CellType | str, value_shape: Sequence[int],
) -> FiniteElement:
    """Create a real element.

    Args:
        cell: Cell to create the element on.
        value_shape: Value shape of the element.

    Returns:
        A 'real' finite element.

    """
    if isinstance(cell, str):
        cell = basix.CellType[cell]

    return RealElement(cell, tuple(value_shape))


def blocked_element(
    sub_element: FiniteElement,
    shape: Sequence[int],
    symmetry: bool | None = None,
) -> FiniteElement:
    """Create a UFLx compatible blocked element.

    Args:
        sub_element: Element used for each block.
        shape: Value shape of the element. For scalar-valued families,
            this can be used to create vector and tensor elements.
        symmetry: Set to ``True`` if the tensor is symmetric. Valid for
            rank 2 elements only.

    Returns:
        A blocked finite element.
    """
    if len(sub_element.reference_value_shape) != 0:
        raise ValueError("Cannot create a blocked element containing a non-scalar element.")

    return BlockedElement(sub_element, shape=shape, symmetry=symmetry)


def wrap_element(element: basix.finite_element.FiniteElement) -> FiniteElement:
    """Wrap a Basix element as a Basix UFLx element."""
    return BasixElement(element)
