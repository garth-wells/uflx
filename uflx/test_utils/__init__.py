# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Utilities for testing UFLx."""

import os
import warnings
from contextlib import contextmanager

from uflx.test_utils.domains import (
    hexahedron,
    interval,
    point,
    quadrilateral,
    tetrahedron,
    triangle,
)
from uflx.test_utils.finite_elements import LagrangeElement

if not os.environ.get("UFLX_ENABLE_TESTING", "0") == "1":
    warnings.warn(
        "Functions in uflx.test_utils are intended to only be used in testing. "
        "Please run with the environment variable `UFLX_ENABLE_TESTING` set to 1 "
        "to disable this warning."
    )
