# Copyright (c) 2024 Matthew Scroggs
# FEniCS Project
# SPDX-License-Identifier: MIT

import numpy as np
import pytest

import basix
import basix_uflx
import uflx


@pytest.mark.parametrize(
    "inputs",
    [
        ("Lagrange", "triangle", 2),
        ("Lagrange", basix.CellType.triangle, 2),
        (basix.ElementFamily.P, basix.CellType.triangle, 2),
        (basix.ElementFamily.P, "triangle", 2),
    ],
)
def test_finite_element(inputs):
    basix_uflx.element(*inputs)


@pytest.mark.parametrize(
    "inputs",
    [
        ("Lagrange", "triangle", 1),
        ("Lagrange", "triangle", 2),
        ("Lagrange", basix.CellType.triangle, 2),
        (basix.ElementFamily.P, basix.CellType.triangle, 2),
        (basix.ElementFamily.P, "triangle", 2),
    ],
)
def test_vector_element(inputs):
    e = basix_uflx.element(*inputs, shape=(2,))
    table = e.tabulate(0, np.array([[0, 0]]))
    assert table.shape == (1, 1, e.dim, e.reference_value_size)


@pytest.mark.parametrize(
    "inputs",
    [
        ("Lagrange", "triangle", 2),
        ("Lagrange", basix.CellType.triangle, 2),
        (basix.ElementFamily.P, basix.CellType.triangle, 2),
        (basix.ElementFamily.P, "triangle", 2),
    ],
)
def test_element(inputs):
    basix_uflx.element(*inputs, shape=(2, 2))


@pytest.mark.parametrize(
    "inputs",
    [
        ("Lagrange", "triangle", 2),
        ("Lagrange", basix.CellType.triangle, 2),
        (basix.ElementFamily.P, basix.CellType.triangle, 2),
        (basix.ElementFamily.P, "triangle", 2),
    ],
)
def test_tensor_element_hash(inputs):
    e = basix_uflx.element(*inputs)
    sym = basix_uflx.blocked_element(e, shape=(2, 2), symmetry=True)
    asym = basix_uflx.blocked_element(e, shape=(2, 2), symmetry=False)
    table = e.tabulate(0, np.array([[0, 0]], dtype=np.float64))
    assert table.shape == (1, 1, e.dim, 1)
    assert sym != asym
    assert hash(sym) != hash(asym)


@pytest.mark.parametrize(
    "cell",
    [
        basix.CellType.triangle,
        basix.CellType.quadrilateral,
        basix.CellType.tetrahedron,
        basix.CellType.prism,
    ],
)
@pytest.mark.parametrize("degree", [1, 3, 6])
@pytest.mark.parametrize("shape", [(), (1,), (2,), (3,), (5,), (2, 2), (3, 3), (4, 1), (5, 1, 7)])
def test_quadrature_element(cell, degree, shape):
    scalar_e = basix_uflx.quadrature_element(cell, (), degree=degree)
    e = basix_uflx.quadrature_element(cell, shape, degree=degree)

    size = 1
    for i in shape:
        size *= i

    assert e.reference_value_size == scalar_e.reference_value_size * size
    assert e.dim == scalar_e.dim * size


@pytest.mark.parametrize(
    "family,cell,degree,shape",
    [
        ("Lagrange", "triangle", 1, None),
        ("Discontinuous Lagrange", "triangle", 1, None),
        ("Lagrange", "quadrilateral", 1, None),
        ("Lagrange", "triangle", 2, None),
        ("Lagrange", "triangle", 1, (2,)),
        ("Lagrange", "triangle", 1, None),
    ],
)
def test_finite_element_eq_hash(family, cell, degree, shape):
    e1 = basix_uflx.element("Lagrange", "triangle", 1, shape=None)
    e2 = basix_uflx.element(family, cell, degree, shape=shape)
    assert (e1 == e2) == (hash(e1) == hash(e2))


@pytest.mark.parametrize(
    "e1,e2",
    [
        (
            basix_uflx.element("Lagrange", "triangle", 1),
            basix_uflx.element("Lagrange", "triangle", 1, shape=(2, 2), symmetry=True),
        ),
        (
            basix_uflx.element("Lagrange", "triangle", 1),
            basix_uflx.element("Lagrange", "triangle", 1, shape=(2, 2)),
        ),
        (
            basix_uflx.element("Lagrange", "triangle", 1),
            basix_uflx.element("Lagrange", "triangle", 1, shape=(2, 2), symmetry=True),
        ),
    ],
)
def test_mixed_element_eq_hash(e1, e2):
    mixed1 = basix_uflx.mixed_element(
        [
            basix_uflx.element("Lagrange", "triangle", 1),
            basix_uflx.element("Lagrange", "triangle", 1, shape=(2, 2), symmetry=True),
        ],
    )
    mixed2 = basix_uflx.mixed_element([e1, e2])
    assert (mixed1 == mixed2) == (hash(mixed1) == hash(mixed2))


@pytest.mark.parametrize(
    ("cell_type", "degree", "reference_map"),
    [
        ("triangle", 2, uflx.maps.IdentityReferenceMap()),
        ("quadrilateral", 2, uflx.maps.IdentityReferenceMap()),
        ("triangle", 3, uflx.maps.IdentityReferenceMap()),
#        ("triangle", 2, basix_uflx._ufl.covariant_piola),
    ],
)
def test_quadrature_element_eq_hash(cell_type, degree, reference_map):
    e1 = basix_uflx.quadrature_element(
        "triangle", scheme="default", degree=2, reference_map=uflx.maps.IdentityReferenceMap(),
    )
    e2 = basix_uflx.quadrature_element(cell_type, scheme="default", degree=degree, reference_map=reference_map)
    assert (e1 == e2) == (hash(e1) == hash(e2))


@pytest.mark.parametrize(
    "cell_type,value_shape", [("triangle", ()), ("quadrilateral", ()), ("triangle", (2,))]
)
def test_real_element_eq_hash(cell_type, value_shape):
    e1 = basix_uflx.real_element("triangle", ())
    e2 = basix_uflx.real_element(cell_type, value_shape)
    assert (e1 == e2) == (hash(e1) == hash(e2))


def test_wrap_element():
    e = basix.create_element(basix.ElementFamily.P, basix.CellType.triangle, 1)
    basix_uflx.wrap_element(e)


def test_dof_ordering():
    e = basix_uflx.element(basix.ElementFamily.P, basix.CellType.triangle, 1)

    e_reordered = basix_uflx.element(
        basix.ElementFamily.P, basix.CellType.triangle, 1, dof_ordering=[1, 2, 0]
    )

    points = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [0.5, 0.5]])

    table = e.tabulate(0, points)[0]
    table2 = e_reordered.tabulate(0, points)[0]

    for i, j in enumerate([1, 2, 0]):
        assert np.allclose(table[:, i], table2[:, j])
