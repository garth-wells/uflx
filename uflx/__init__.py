# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""UFLx: Unified Form Language."""

from uflx.domains import domain
from uflx.function_spaces import function_space
from uflx.functions import TestFunction, TrialFunction
from uflx.integrals import dS, ds, dx
from uflx.operators import inner, grad
