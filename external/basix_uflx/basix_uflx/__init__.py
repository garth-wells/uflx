"""Basix UFLx interface."""

from abc import abstractmethod
from collections.abc import Sequence

import hashlib as _hashlib
import itertools as _itertools
from abc import abstractmethod as _abstractmethod
from warnings import warn as _warn

import numpy as np
import numpy.typing as _npt
import uflx as _uflx
from uflx.finite_elements import AbstractReferenceMappedFiniteElement as _ARMFE
from uflx.scalars import AbstractInteger as _AbstractInteger

import basix as _basix

__all__ = [
    "element",
    "enriched_element",
    "custom_element",
    "mixed_element",
    "quadrature_element",
    "real_element",
    "blocked_element",
    "wrap_element",
]


def convert_map(basix_map: _basix.MapType) -> _uflx.maps.AbstractReferenceMap:
    """Convert a basix map to a UFLx map."""
    match basix_map:
        case _basix.MapType.identity:
            return _uflx.maps.IdentityReferenceMap()
        case _:
            raise NotImplementedError()


class BasixCell(_uflx.entities.AbstractEntity):
    """A cell defined by Basix."""

    def __init__(self, basix_cell: _basix.CellType):
        self._basix_cell = basix_cell

    def __eq__(self, other) -> bool:
        """Check if this entity is equal to another entity."""
        return isinstance(other, BasixCell) and self._basix_cell == other.basix_cell

    @property
    def topological_dimension(self) -> int:
        """Topological dimension of the entity."""
        return len(_basix.cell.topology(self._basix_cell)) - 1

    def sub_entities(self, dim: int) -> list[_uflx.entities.AbstractEntity]:
        """Get a list of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            A list of sub-entities of the given dimension.
        """
        return _basix.cell.subentity_types(self._basix_cell)[dim]

    def __hash__(self):
        """Hash."""
        return hash(("basix.uflx", self._basix_cell))


# TODO: move this base class into uflx_codegeneration
class _ElementBase(_ARMFE):
    """A base wrapper to allow elements to be used with UFLx.

    This class includes methods and properties needed by UFLx and FFCx.
    This is a base class containing functions common to all the element
    types defined in this file.
    """

    def __init__(self, repr: str):
        """Initialise the element."""
        self._repr = repr

    def __hash__(self):
        """Hash."""
        return hash(("basix.uflx", self._repr))

    def __repr__(self) -> str:
        """Representation."""
        return self._repr

    def __str__(self):
        """String representation."""
        return self._repr

    @property
    @abstractmethod
    def dim(self) -> int | _AbstractInteger:
        """The dimension of the finite element, ie the number of basis functions."""

    @abstractmethod
    def tabulate(self, points: np.ndarray, derivative: tuple[int, ...]) -> np.ndarray:
        """Create table of basis function values.

        TODO: document return shape
        """


class _BasixElement(_ElementBase):
    """A wrapper allowing Basix elements to be used directly with UFLx.

    This class allows elements created with `basix.create_element` to be
    wrapped as UFLx compatible elements. Users should not directly call
    this class's initialiser, but should use the `element` function
    instead.
    """

    _element: _basix.finite_element.FiniteElement

    def __init__(self, element: _basix.finite_element.FiniteElement):
        """Create a Basix element."""
        if element.family == _basix.ElementFamily.custom:
            repr = f"custom Basix element ({_compute_signature(element)})"
        else:
            repr = (
                f"Basix element ({element.family.name}, {element.cell_type.name}, "
                f"{element.degree}, "
                f"{element.lagrange_variant.name}, {element.dpc_variant.name}, "
                f"{element.discontinuous}, "
                f"{element.dtype}, {element.dof_ordering})"
            )

        super().__init__(repr)
        self._element = element

    def __eq__(self, other) -> bool:
        return isinstance(other, _BasixElement) and self._element == other._element

    @property
    def cell(self) -> _uflx.entities.AbstractEntity:
        return BasixCell(self._element.cell_type)

    @property
    def dim(self) -> int:
        return self._element.dim

    @property
    def lagrange_superdegree(self) -> int:
        return self._element.embedded_superdegree

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        raise NotImplementedError()

    @property
    def reference_map(self) -> _uflx.maps.AbstractReferenceMap:
        return convert_map(self._element.map_type)

    @property
    def reference_value_shape(self) -> tuple[int, ...]:
        return self._element.value_shape

    def tabulate(self, points: np.ndarray, derivative: tuple[int, ...]) -> np.ndarray:
        return self._element.tabulate(sum(derivative), points)[derivative]


class _ComponentElement(_ElementBase):
    """An element representing one component of a _BasixElement.

    This element type is used when UFLx's ``get_component_element``
    function is called.

    """

    _element: _ElementBase
    _component: int

    def __init__(self, element: _ElementBase, component: int):
        """Initialise the element."""
        self._element = element
        self._component = component
        repr = f"component element ({element!r}, {component}"
        repr += ")"
        super().__init__(repr, element.cell_type.name, (1,), element._degree)


class _MixedElement(_ElementBase):
    """A mixed element that combines two or more elements.

    This can be used when multiple different elements appear in a form.
    Users should not directly call this class's initilizer, but should
    use the :func:`mixed_element` function instead.
    """

    _sub_elements: list[_ElementBase]

    def __init__(self, sub_elements: list[_ElementBase]):
        """Initialise the element."""
        assert len(sub_elements) > 0
        self._sub_elements = sub_elements
        reference_map = (
            _uflx.identity_pullback
            if all(isinstance(e.pullback, _IdentityPullback) for e in sub_elements)
            else _MixedPullback(self)
        )

        super().__init__("mixed element (" + ", ".join(i._repr for i in sub_elements) + ")")


class _BlockedElement(_ElementBase):
    """Element with a block size that contains multiple copies of a sub element.

    This can be used to (for example) create vector and tensor Lagrange
    elements. Users should not directly call this classes initilizer,
    but should use the `blocked_element` function instead.

    """

    _block_shape: tuple[int, ...]
    _sub_element: _ElementBase
    _block_size: int
    _has_symmetry: bool

    def __init__(
        self,
        sub_element: _ElementBase,
        shape: tuple[int, ...],
        symmetry: bool | None = None,
    ):
        """Initialise the element."""
        if sub_element.reference_value_size != 1:
            raise ValueError(
                "Blocked elements of non-scalar elements are not supported. "
                "Try using _MixedElement instead."
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

        repr = f"blocked element ({sub_element!r}, {shape}"
        if symmetry is not None:
            repr += f", symmetry={symmetry}"
        repr += ")"

        super().__init__(repr)

        if symmetry:
            n = 0
            symmetry_mapping: dict[tuple[int, ...], int] = {}
            for i in range(shape[0]):
                for j in range(i + 1):
                    symmetry_mapping[(i, j)] = n
                    symmetry_mapping[(j, i)] = n
                    n += 1

            self._reference_map = _SymmetricPullback(self, symmetry_mapping)


class _QuadratureElement(_ElementBase):
    """A quadrature element."""

    def __init__(
        self,
        cell: _basix.CellType,
        points: _npt.NDArray[np.floating],
        weights: _npt.NDArray[np.floating],
        reference_map: _uflx.maps.AbstractReferenceMap,
        degree: int | None = None,
        dtype: _npt.DTypeLike = np.float64,
    ):
        """Initialise the element."""
        self._points = points.astype(dtype)
        self._weights = weights.astype(dtype)
        repr = f"QuadratureElement({cell.name}, {points!r}, {weights!r}, {reference_map!r})".replace(
            "\n", ""
        )
        self._cell_type = cell
        self._entity_counts = [len(i) for i in _basix.topology(cell)]
        self._reference_map = reference_map

        if degree is None:
            degree = len(points)

        super().__init__(repr)


class _RealElement(_ElementBase):
    """A real element."""

    def __init__(self, cell: _basix.CellType, value_shape: tuple[int, ...]):
        """Initialise the element."""
        self._cell_type = cell
        tdim = len(_basix.topology(cell)) - 1

        super().__init__(f"RealElement({cell.name}, {value_shape})", cell.name, value_shape, 0)

        self._entity_counts = []
        if tdim >= 1:
            self._entity_counts.append(self.cell.num_vertices)
        if tdim >= 2:
            self._entity_counts.append(self.cell.num_edges)
        if tdim >= 3:
            self._entity_counts.append(self.cell.num_facets)
        self._entity_counts.append(1)


def element(
    family: _basix.ElementFamily | str,
    cell: _basix.CellType | str,
    degree: int,
    lagrange_variant: _basix.LagrangeVariant = _basix.LagrangeVariant.unset,
    dpc_variant: _basix.DPCVariant = _basix.DPCVariant.unset,
    discontinuous: bool = False,
    shape: tuple[int, ...] | None = None,
    symmetry: bool | None = None,
    dof_ordering: list[int] | None = None,
    dtype: _npt.DTypeLike | None = None,
) -> _ElementBase:
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
        cell = _basix.CellType[cell]
    if isinstance(family, str):
        if family.startswith("Discontinuous "):
            family = family[14:]
            discontinuous = True
        if family in ["DP", "DG", "DQ"]:
            family = "P"
            discontinuous = True
        if family == "CG":
            _warn(
                '"CG" element name is deprecated. Consider using "Lagrange" or "P" instead',
                DeprecationWarning,
                stacklevel=2,
            )
            family = "P"
            discontinuous = False
        if family == "DPC":
            discontinuous = True

        family = _basix.finite_element.string_to_family(family, cell.name)

    # Default variant choices
    EF = _basix.ElementFamily
    if lagrange_variant == _basix.LagrangeVariant.unset:
        if family == EF.P:
            lagrange_variant = _basix.LagrangeVariant.gll_warped
        elif family in [EF.RT, EF.N1E]:
            lagrange_variant = _basix.LagrangeVariant.legendre
        elif family in [EF.serendipity, EF.BDM, EF.N2E]:
            lagrange_variant = _basix.LagrangeVariant.legendre

    if dpc_variant == _basix.DPCVariant.unset:
        if family in [EF.serendipity, EF.BDM, EF.N2E]:
            dpc_variant = _basix.DPCVariant.legendre
        elif family == EF.DPC:
            dpc_variant = _basix.DPCVariant.diagonal_gll

    e = _basix.create_element(
        family,
        cell,
        degree,
        lagrange_variant,
        dpc_variant,
        discontinuous,
        dof_ordering=dof_ordering,
        dtype=dtype,
    )
    ufl_e = _BasixElement(e)

    if shape is None or shape == tuple(e.value_shape):
        if symmetry is not None:
            raise ValueError("Cannot pass a symmetry argument to this element.")
        return ufl_e
    else:
        return blocked_element(ufl_e, shape=shape, symmetry=symmetry)


def enriched_element(
    elements: list[_ElementBase],
    map_type: _basix.MapType | None = None,
) -> _ElementBase:
    """Create an UFLx compatible enriched element from a list of elements.

    Args:
        elements: The list of elements
        map_type: The map type for the enriched element.

    Returns:
        An enriched finite element.

    """
    ct = elements[0].cell_type
    ptype = elements[0].polyset_type
    vshape = elements[0].reference_value_shape
    vsize = elements[0].reference_value_size
    if map_type is None:
        map_type = elements[0].map_type
        for e in elements:
            if e.map_type != map_type:
                raise ValueError("Enriched elements on different map types not supported.")

    dtype = e.dtype
    hcd = min(e.embedded_subdegree for e in elements)
    hd = max(e.embedded_superdegree for e in elements)
    ss = _basix.sobolev_spaces.intersection([e.basix_sobolev_space for e in elements])
    discontinuous = True
    for e in elements:
        if not e.discontinuous:
            discontinuous = False
        if e.cell_type != ct:
            raise ValueError("Enriched elements on different cell types not supported.")
        if e.polyset_type != ptype:
            raise ValueError("Enriched elements on different polyset types not supported.")
        if e.reference_value_shape != vshape or e.reference_value_size != vsize:
            raise ValueError("Enriched elements on different value shapes not supported.")
        if e.dtype != dtype:
            raise ValueError("Enriched elements with different dtypes no supported.")
    nderivs = max(e.interpolation_nderivs for e in elements)

    x = []
    for pts_lists in zip(*[e._x for e in elements]):
        x.append([np.concatenate(pts) for pts in zip(*pts_lists)])
    M = []
    for M_lists in zip(*[e._M for e in elements]):
        M_row = []
        for M_parts in zip(*M_lists):
            ndofs = sum(mat.shape[0] for mat in M_parts)
            npts = sum(mat.shape[2] for mat in M_parts)
            deriv_dim = max(mat.shape[3] for mat in M_parts)
            new_M = np.zeros((ndofs, vsize, npts, deriv_dim))
            pt = 0
            dof = 0
            for mat in M_parts:
                new_M[dof : dof + mat.shape[0], :, pt : pt + mat.shape[2], : mat.shape[3]] = mat
                dof += mat.shape[0]
                pt += mat.shape[2]
            M_row.append(new_M)
        M.append(M_row)

    dim = sum(e.dim for e in elements)
    wcoeffs = np.zeros(
        (dim, _basix.polynomials.dim(_basix.PolynomialType.legendre, ct, hd) * vsize)
    )
    row = 0
    for e in elements:
        wcoeffs[row : row + e.dim, :] = _basix.polynomials.reshape_coefficients(
            _basix.PolynomialType.legendre,
            ct,
            e._wcoeffs,  # type: ignore
            vsize,
            e.embedded_superdegree,
            hd,
        )
        row += e.dim

    return custom_element(
        ct,
        list(vshape),
        wcoeffs,
        x,
        M,
        nderivs,
        map_type,
        ss,
        discontinuous,
        hcd,
        hd,
        ptype,
        dtype=dtype,
    )


def custom_element(
    cell_type: _basix.CellType,
    reference_value_shape: Sequence[int],
    wcoeffs: _npt.ArrayLike,
    x: Sequence[Sequence[_npt.ArrayLike]],
    M: Sequence[Sequence[_npt.ArrayLike]],
    interpolation_nderivs: int,
    map_type: _basix.MapType,
    sobolev_space: _basix.SobolevSpace,
    discontinuous: bool,
    embedded_subdegree: int,
    embedded_superdegree: int,
    polyset_type: _basix.PolysetType = _basix.PolysetType.standard,
    dtype: _npt.DTypeLike | None = None,
) -> _ElementBase:
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
    e = _basix.create_custom_element(
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
    return _BasixElement(e)


def mixed_element(elements: Sequence[_ElementBase]) -> _ElementBase:
    """Create a UFLx compatible mixed element from a list of elements.

    Args:
        elements: The list of elements

    Returns:
        A mixed finite element.
    """
    return _MixedElement(list(elements))


def quadrature_element(
    cell: str | _basix.CellType,
    value_shape: Sequence[int] = (),
    scheme: str | None = None,
    degree: int | None = None,
    points: _npt.NDArray[np.floating] | None = None,
    weights: _npt.NDArray[np.floating] | None = None,
    reference_map: _uflx.maps.AbstractReferenceMap | None = None,
    symmetry: bool | None = None,
    dtype: _npt.DTypeLike | None = None,
) -> _ElementBase:
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
        cell = _basix.CellType[cell]

    if reference_map is None:
        refernece_map = _uflx.maps.IdentityReferenceMap()

    if points is None:
        assert weights is None
        assert degree is not None
        if scheme is None:
            points, weights = _basix.make_quadrature(cell, degree)  # type: ignore
        else:
            points, weights = _basix.make_quadrature(  # type: ignore
                cell, degree, rule=_basix.quadrature.string_to_type(scheme)
            )

    assert points is not None
    assert weights is not None

    e = _QuadratureElement(cell, points, weights, reference_map, degree, dtype=dtype)
    if tuple(value_shape) == ():
        if symmetry is not None:
            raise ValueError("Cannot pass a symmetry argument to this element.")
        return e
    else:
        return blocked_element(e, shape=tuple(value_shape), symmetry=symmetry)


def real_element(
    cell: _basix.CellType | str, value_shape: Sequence[int],
) -> _ElementBase:
    """Create a real element.

    Args:
        cell: Cell to create the element on.
        value_shape: Value shape of the element.

    Returns:
        A 'real' finite element.

    """
    if isinstance(cell, str):
        cell = _basix.CellType[cell]

    return _RealElement(cell, tuple(value_shape))


def blocked_element(
    sub_element: _ElementBase,
    shape: Sequence[int],
    symmetry: bool | None = None,
) -> _ElementBase:
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

    return _BlockedElement(sub_element, shape=shape, symmetry=symmetry)


def wrap_element(element: _basix.finite_element.FiniteElement) -> _ElementBase:
    """Wrap a Basix element as a Basix UFLx element."""
    return _BasixElement(element)
