# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test entities."""


def test_euler_characteristic(entity):
    """Test Euler characteristic of entity."""
    match entity.topological_dimension:
        case 0:
            assert len(entity.sub_entities(0)) == 1
        case 1:
            assert len(entity.sub_entities(0)) - len(entity.sub_entities(1)) == 1
        case 2:
            assert (
                len(entity.sub_entities(0))
                - len(entity.sub_entities(1))
                + len(entity.sub_entities(2))
                == 1
            )
        case 3:
            assert (
                len(entity.sub_entities(0))
                - len(entity.sub_entities(1))
                + len(entity.sub_entities(2))
                == 2
            )


def test_is_simplex(entity):
    """Test that is_simplex correctly identifies simplex vs tensor-product cells."""
    simplex_names = {"point", "interval", "triangle", "tetrahedron"}
    assert entity.is_simplex == (entity.name in simplex_names)
