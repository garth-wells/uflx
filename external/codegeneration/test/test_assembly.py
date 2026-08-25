"""Test code generation."""

import os

import numpy as np
import pytest
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

code_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), ".code")
if not os.path.isdir(code_dir):
    os.mkdir(code_dir)


@pytest.mark.parametrize(
    ("cell", "expected_mat"),
    [
        pytest.param(
            [0, 1, 3],
            np.array(
                [
                    [0.025, 0.0125, 0.0125],
                    [0.0125, 0.025, 0.0125],
                    [0.0125, 0.0125, 0.025],
                ]
            ),
            id="cell_0",
        ),
        pytest.param(
            [1, 4, 3],
            np.array(
                [
                    [0.025, 0.0125, 0.0125],
                    [0.0125, 0.025, 0.0125],
                    [0.0125, 0.0125, 0.025],
                ]
            ),
            id="cell_1",
        ),
        pytest.param(
            [1, 2, 4],
            np.array(
                [
                    [7 / 120, 7 / 240, 7 / 240],
                    [7 / 240, 7 / 120, 7 / 240],
                    [7 / 240, 7 / 240, 7 / 120],
                ]
            ),
            id="cell_2",
        ),
        pytest.param(
            [2, 5, 4],
            np.array(
                [
                    [7 / 120, 7 / 240, 7 / 240],
                    [7 / 240, 7 / 120, 7 / 240],
                    [7 / 240, 7 / 240, 7 / 120],
                ]
            ),
            id="cell_3",
        ),
    ],
)
def test_mass_matrix(lagrange_element, cell, expected_mat):
    """Test code generation for a mass matrix."""
    element = lagrange_element("triangle", 1)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    code, signatures = uflx_codegeneration.generate(form)

    pts = np.array([[0.0, 0.0], [0.3, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 1.0], [1.0, 1.0]])

    filename = "test_mass_matrix"

    ffi = FFI()
    ffi.cdef("\n".join(signatures.values()))
    ffi.set_source(f"{filename}", code)
    so = ffi.compile(code_dir)

    lib = ffi.dlopen(so)

    mat = np.zeros((3, 3))
    coords = np.zeros((3, 2))
    empty = np.zeros(0)

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


@pytest.mark.parametrize(
    ("cell", "expected_mat"),
    [
        pytest.param(
            [0, 1, 3],
            np.array(
                [
                    [1.81666666666667, -1.66666666666667, -0.15],
                    [-1.66666666666667, 1.66666666666667, 0],
                    [-0.15, 0, 0.15],
                ]
            ),
            id="cell_0",
        ),
        pytest.param(
            [1, 4, 3],
            np.array(
                [
                    [0.15, -0.15, 0],
                    [-0.15, 1.81666666666667, -1.66666666666667],
                    [0, -1.66666666666667, 1.66666666666667],
                ]
            ),
            id="cell_1",
        ),
        pytest.param(
            [1, 2, 4],
            np.array(
                [
                    [1.06428571428571, -0.714285714285714, -0.35],
                    [-0.714285714285714, 0.714285714285714, 0],
                    [-0.35, 0, 0.35],
                ]
            ),
            id="cell_2",
        ),
        pytest.param(
            [2, 5, 4],
            np.array(
                [
                    [0.35, -0.35, 0],
                    [-0.35, 1.06428571428571, -0.714285714285714],
                    [0, -0.714285714285714, 0.714285714285714],
                ]
            ),
            id="cell_3",
        ),
    ],
)
def test_stiffness_matrix(lagrange_element, cell, expected_mat):
    """Test code generation for a stiffness matrix."""
    element = lagrange_element("triangle", 1)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(grad(u), grad(v)) * dx

    code, signatures = uflx_codegeneration.generate(form)

    pts = np.array([[0.0, 0.0], [0.3, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 1.0], [1.0, 1.0]])

    filename = "test_stiffness_matrix"

    ffi = FFI()
    ffi.cdef("\n".join(signatures.values()))
    ffi.set_source(f"{filename}", code)
    so = ffi.compile(code_dir)

    lib = ffi.dlopen(so)

    mat = np.zeros((3, 3))
    coords = np.zeros((3, 2))
    empty = np.zeros(0)

    for i, j in enumerate(cell):
        for k, p in enumerate(pts[j]):
            coords[i, k] = p
    print(cell)
    print(coords)
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
    print(mat)
    print(expected_mat)
    assert np.allclose(mat, expected_mat)


@pytest.mark.parametrize(
    ("cell", "expected_vec"),
    [
        pytest.param(
            [0, 1, 3],
            np.array([0.0037500000000000033, 0.0075, 0.0037500000000000033]),
            id="cell_0",
        ),
        pytest.param(
            [1, 4, 3],
            np.array([0.011249999999999982, 0.01125, 0.0075]),
            id="cell_1",
        ),
        pytest.param(
            [1, 2, 4],
            np.array([0.05541666666666667, 0.07583333333333332, 0.05541666666666664]),
            id="cell_2",
        ),
        pytest.param(
            [2, 5, 4],
            np.array([0.09625000000000011, 0.09624999999999997, 0.07583333333333334]),
            id="cell_3",
        ),
    ],
)
def test_linear_form(lagrange_element, cell, expected_vec):
    """Test code generation for a mass matrix."""
    element = lagrange_element("triangle", 1)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    v = TestFunction(space)
    x = SpatialCoordinate(2)
    form = x[0] * v * dx

    code, signatures = uflx_codegeneration.generate(form)

    pts = np.array([[0.0, 0.0], [0.3, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 1.0], [1.0, 1.0]])

    filename = "test_linear_form"

    ffi = FFI()
    ffi.cdef("\n".join(signatures.values()))
    ffi.set_source(f"{filename}", code)
    so = ffi.compile(code_dir)

    lib = ffi.dlopen(so)

    vec = np.zeros(3)
    coords = np.zeros((3, 2))
    empty = np.zeros(0)

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
