# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test finite elements."""

import pytest
from utils import LagrangeElement, hexahedron, interval, point, quadrilateral, tetrahedron, triangle


@pytest.mark.parametrize(
    "cell",
    [
        point,
        interval,
        triangle,
        quadrilateral,
        tetrahedron,
        hexahedron,
    ],
)
def test_lagrange_element(cell):
    """Test Lagrange element properties."""
    element = LagrangeElement(cell, 2)

    assert element.cell == cell
    assert element.reference_value_shape == ()
    assert element.embedded_lagrange_superdegree == 2
