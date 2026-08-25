"""Test code generation."""

import numpy as np
from cffi import FFI
from uflx import (
    SpatialCoordinate,
    TestFunction,
    TrialFunction,
    coordinate_element,
    dx,
    function_space,
    grad,
    inner,
)

import uflx_codegeneration


def test_mass_matrix(lagrange_element, code_dir):
    """Test code generation for a mass matrix."""
    element = lagrange_element("triangle", 1)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    code, signatures = uflx_codegeneration.generate(form)

    pts = np.array([[0.0, 0.0], [0.3, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 1.0], [1.0, 1.0]])

    cells = [[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]]
    expected_matrices = [
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
    coords = np.zeros((3, 2))
    empty = np.zeros(0)

    for cell, expected_mat in zip(cells, expected_matrices):
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


def test_stiffness_matrix(lagrange_element, code_dir):
    """Test code generation for a stiffness matrix."""
    element = lagrange_element("triangle", 1)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(grad(u), grad(v)) * dx

    code, signatures = uflx_codegeneration.generate(form)

    pts = np.array([[0.0, 0.0], [0.3, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 1.0], [1.0, 1.0]])
    cells = [[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]]
    expected_matrices = [
        np.array(
            [
                [1.81666666666667, -1.66666666666667, -0.15],
                [-1.66666666666667, 1.66666666666667, 0],
                [-0.15, 0, 0.15],
            ]
        ),
        np.array(
            [
                [0.15, -0.15, 0],
                [-0.15, 1.81666666666667, -1.66666666666667],
                [0, -1.66666666666667, 1.66666666666667],
            ]
        ),
        np.array(
            [
                [1.06428571428571, -0.714285714285714, -0.35],
                [-0.714285714285714, 0.714285714285714, 0],
                [-0.35, 0, 0.35],
            ]
        ),
        np.array(
            [
                [0.35, -0.35, 0],
                [-0.35, 1.06428571428571, -0.714285714285714],
                [0, -0.714285714285714, 0.714285714285714],
            ]
        ),
    ]

    filename = "test_stiffness_matrix"

    ffi = FFI()
    ffi.cdef("\n".join(signatures.values()))
    ffi.set_source(f"{filename}", code)
    so = ffi.compile(code_dir)

    lib = ffi.dlopen(so)

    mat = np.zeros((3, 3))
    coords = np.zeros((3, 2))
    empty = np.zeros(0)

    for cell, expected_mat in zip(cells, expected_matrices):
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


def test_linear_form(lagrange_element, code_dir):
    """Test code generation for a mass matrix."""
    element = lagrange_element("triangle", 1)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    v = TestFunction(space)
    x = SpatialCoordinate(2)
    form = x[0] * v * dx

    code, signatures = uflx_codegeneration.generate(form)

    pts = np.array([[0.0, 0.0], [0.3, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 1.0], [1.0, 1.0]])
    cells = [[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]]

    expected_vectors = [
        np.array([0.0037500000000000033, 0.0075, 0.0037500000000000033]),
        np.array([0.011249999999999982, 0.01125, 0.0075]),
        np.array([0.05541666666666667, 0.07583333333333332, 0.05541666666666664]),
        np.array([0.09625000000000011, 0.09624999999999997, 0.07583333333333334]),
    ]

    filename = "test_linear_form"

    ffi = FFI()
    ffi.cdef("\n".join(signatures.values()))
    ffi.set_source(f"{filename}", code)
    so = ffi.compile(code_dir)

    lib = ffi.dlopen(so)

    vec = np.zeros(3)
    coords = np.zeros((3, 2))
    empty = np.zeros(0)

    for cell, expected_vec in zip(cells, expected_vectors):
        for i, j in enumerate(cell):
            for k, p in enumerate(pts[j]):
                coords[i, k] = p
        vec[:] = 0.0
        lib.tabulate_tensor_f64(
            ffi.cast("double*", vec.ctypes.data),
            ffi.cast("double*", empty.ctypes.data),
            ffi.cast("double*", empty.ctypes.data),
            ffi.cast("double*", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        assert np.allclose(vec, expected_vec)
