# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test function spaces."""

import pytest

from uflx import coordinate_element, function_space


@pytest.mark.parametrize(
    ("cell_name", "gdim"),
    [
        ("triangle", 2),
        ("triangle", 3),
        ("quadrilateral", 2),
        ("quadrilateral", 3),
        ("tetrahedron", 3),
    ],
)
def test_function_space(cell_name, gdim, lagrange_element, triangle, quadrilateral, tetrahedron):
    """Test function space with single cell."""
    cell = {"triangle": triangle, "quadrilateral": quadrilateral, "tetrahedron": tetrahedron}[cell_name]
    domain = coordinate_element(lagrange_element(cell, 1, (3,)))
    space = function_space(domain, lagrange_element(cell, 1))
    assert len(space.elements) == len(space.domain.cells) == 1


@pytest.mark.parametrize("gdim", [2, 3])
def test_function_space_multiple_cells(gdim, lagrange_element, triangle, quadrilateral):
    """Test function space with multiple cells."""
    domain = coordinate_element(
        [
            lagrange_element(triangle, 1, (3,)),
            lagrange_element(quadrilateral, 1, (3,)),
        ]
    )
    space = function_space(
        domain,
        [lagrange_element(triangle, 2), lagrange_element(quadrilateral, 2)],
    )
    assert len(space.elements) == len(space.domain.cells) == 2
