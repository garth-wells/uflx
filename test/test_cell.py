# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test cells."""

import pytest
from utils import interval, point, tetrahedron, triangle


@pytest.mark.parametrize(
    "cell",
    [
        point,
        interval,
        triangle,
        tetrahedron,
    ],
)
def test_euler_characteristic(cell):
    """Test Euler characteristic of cell."""
    match cell.topological_dimension:
        case 0:
            assert len(cell.sub_entities(0)) == 1
        case 1:
            assert len(cell.sub_entities(0)) - len(cell.sub_entities(1)) == 1
        case 2:
            assert (
                len(cell.sub_entities(0)) - len(cell.sub_entities(1)) + len(cell.sub_entities(2))
                == 1
            )
        case 3:
            assert (
                len(cell.sub_entities(0)) - len(cell.sub_entities(1)) + len(cell.sub_entities(2))
                == 2
            )
