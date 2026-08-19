# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test finite elements."""

import pytest


def test_lagrange_element(entity, lagrange_element):
    """Test Lagrange element properties."""
    element = lagrange_element(entity, 2)

    assert element.cell == entity
    assert element.reference_value_shape == ()
    assert element.lagrange_superdegree == 2
