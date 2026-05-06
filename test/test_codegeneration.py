"""Test code generation."""

import numpy as np
from utils import LagrangeElement, triangle

from uflx import TestFunction, TrialFunction, domain, dx, function_space, inner
from uflx import codegeneration


def test_mass_matrix():
    """Test code generation for a mass matrix."""
    element = LagrangeElement(triangle, 2)
    space = function_space(domain(triangle), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    code = codegeneration.generate(form)

    print(code)

    pts = np.array([[0., 0.], [0.3, 0.], [1., 0.], [0., 1.], [0.3, 1.], [1., 1.]])
    cells = np.array([[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]])

    expected_local_matrices = [
        np.array([
            [0.025, 0.0125000000000000 0.0125000000000000],
            [0.0125, 0.0250000000000000 0.0125000000000000],
            [0.0125, 0.0125000000000000 0.0250000000000000]
        ]),
        np.array([
local matrix for [1, 4, 3]
0.0250000000000000 0.0125000000000000 0.0125000000000000 
0.0125000000000000 0.0250000000000000 0.0125000000000000 
0.0125000000000000 0.0125000000000000 0.0250000000000000 
        ]),
        np.array([

local matrix for [1, 2, 4]
0.0583333333333333 0.0291666666666666 0.0291666666666667 
0.0291666666666666 0.0583333333333333 0.0291666666666666 
0.0291666666666667 0.0291666666666666 0.0583333333333333 
        ]),
        np.array([

local matrix for [2, 5, 4]
0.0583333333333332 0.0291666666666668 0.0291666666666667 
0.0291666666666668 0.0583333333333333 0.0291666666666666 
0.0291666666666667 0.0291666666666666 0.0583333333333333
        ])
    ]
    expected_matrix = np.array([
        [0.025, 0.0125, 0., 0.0125, 0., 0.],
        [0.0125, 13/120, 7/240, 0.025, 1/24, 0.],
        [0., 7/240, 7/60, 0., 7/120, 7/240],
        [0.0125, 0.025, 0., 0.05, 0.0125, 0.],
        [0., 1/24, 7/120, 0.0125, 17/120, 7/240],
        [0., 0., 7/240, 0, 7/240, 7/120],
    ])
