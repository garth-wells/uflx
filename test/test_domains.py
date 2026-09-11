# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test domains."""

from uflx import coordinate_element


def test_is_affine_map_true_for_degree_one_simplex(lagrange_element):
    """A degree 1 Lagrange coordinate element on a simplex cell is an affine map."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    assert domain.is_affine_map

    domain = coordinate_element(lagrange_element("tetrahedron", 1, (3,)))
    assert domain.is_affine_map


def test_is_affine_map_false_for_higher_degree(lagrange_element):
    """A higher degree Lagrange coordinate element is not an affine map."""
    domain = coordinate_element(lagrange_element("triangle", 2, (2,)))
    assert not domain.is_affine_map


def test_is_affine_map_false_for_tensor_product_cell(lagrange_element):
    """A degree 1 Lagrange coordinate element on a non-simplex cell is not affine.

    Even at degree 1, a quadrilateral or hexahedron's coordinate map is
    multilinear, not affine, since these cells aren't simplices.
    """
    domain = coordinate_element(lagrange_element("quadrilateral", 1, (2,)))
    assert not domain.is_affine_map

    domain = coordinate_element(lagrange_element("hexahedron", 1, (3,)))
    assert not domain.is_affine_map
