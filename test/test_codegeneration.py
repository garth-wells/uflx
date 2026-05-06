"""Test code generation."""

import os

import numpy as np
from cffi import FFI
from utils import LagrangeElement, triangle

from uflx import TestFunction, TrialFunction, codegeneration, domain, dx, function_space, inner

code_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), ".code")
if not os.path.isdir(code_dir):
    os.mkdir(code_dir)


def test_mass_matrix():
    """Test code generation for a mass matrix."""
    element = LagrangeElement(triangle, 2)
    space = function_space(domain(triangle), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    code, signatures = codegeneration.generate(form)

    pts = np.array([[0.0, 0.0], [0.3, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 1.0], [1.0, 1.0]])
    cells = np.array([[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]])

    expected_local_matrices = [
        np.array(
            [
                [0.025, 0.0125, 0.0125],
                [0.0125, 0.025, 0.0125],
                [0.0125, 0.0125, 0.025],
            ]
        ),
        np.array(
            [
                [0.025, 0.0125, 0.0125],
                [0.0125, 0.025, 0.0125],
                [0.0125, 0.0125, 0.025],
            ]
        ),
        np.array(
            [
                [7 / 120, 7 / 240, 7 / 240],
                [7 / 240, 7 / 120, 7 / 240],
                [7 / 240, 7 / 240, 7 / 120],
            ]
        ),
        np.array(
            [
                [7 / 120, 7 / 240, 7 / 240],
                [7 / 240, 7 / 120, 7 / 240],
                [7 / 240, 7 / 240, 7 / 120],
            ]
        ),
    ]

    filename = "test_mass_matrix"

    ffi = FFI()
    ffi.cdef("\n".join(signatures.values()))
    ffi.set_source(f"{filename}", code)
    so = ffi.compile(code_dir)

    lib = ffi.dlopen(so)

    mat = np.zeros((3, 3))
    coords = np.zeros((3, 3))
    empty = np.zeros(0)

    for cell, expected_mat in zip(cells, expected_local_matrices):
        for i, j in enumerate(cell):
            for k, p in enumerate(pts[j]):
                coords[i, k] = p
        mat[:] = 0.0
        lib.tabulate_tensor_f64(
            ffi.cast("double*", mat.ctypes.data),
            ffi.cast("double*", empty.ctypes.data),
            ffi.cast("double*", empty.ctypes.data),
            ffi.cast("double*", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        assert np.allclose(mat, expected_mat)
