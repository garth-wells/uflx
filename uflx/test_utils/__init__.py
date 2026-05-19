# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Utilities for testing UFLx."""

from contextlib import contextmanager
import os
import warnings

from uflx.test_utils.domains import interval, triangle, quadrilateral, tetrahedron, hexahedron, point
from uflx.test_utils.finite_elements import LagrangeElement

if not os.environ.get("UFLX_ENABLE_TESTING", "0") == "1":
    warnings.warn(
        "Functions in uflx.test_utils are intended to only be used in testing. "
        "Please run with the environment variable `UFLX_ENABLE_TESTING` set to 1 "
        "to disable this warning."
    )
