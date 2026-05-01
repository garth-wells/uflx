# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test function spaces."""

import pytest
from utils import LagrangeElement, quadrilateral, tetrahedron, triangle

from uflx import domain, function_space


@pytest.mark.parametrize(
    ("cell", "gdim"),
    [
        (triangle, 2),
        (triangle, 3),
        (quadrilateral, 2),
        (quadrilateral, 3),
        (tetrahedron, 3),
    ],
)
def test_function_space(cell, gdim):
    """Test function space with single cell."""
    space = function_space(domain(cell, gdim), LagrangeElement(cell, 1))
    assert len(space.elements) == len(space.domain.cells) == 1


@pytest.mark.parametrize("gdim", [2, 3])
def test_function_space_multiple_cells(gdim):
    """Test function space with multiple cells."""
    space = function_space(
        domain([triangle, quadrilateral], gdim),
        [LagrangeElement(triangle, 2), LagrangeElement(quadrilateral, 2)],
    )
    assert len(space.elements) == len(space.domain.cells) == 2
