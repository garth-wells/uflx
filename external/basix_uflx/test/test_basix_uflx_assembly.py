"""Test code generation."""

import os

import numpy as np
import pytest
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

from basix_uflx import element

code_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), ".code")
if not os.path.isdir(code_dir):
    os.mkdir(code_dir)


@pytest.mark.parametrize(
    ("element_degree", "geometry_degree", "coords", "expected_local_matrix"),
    [
        pytest.param(
            1,
            1,
            np.array([[0.0, 0.0], [0.3, 0.0], [0.0, 1.0]]),
            np.array(
                [
                    [0.025, 0.0125, 0.0125],
                    [0.0125, 0.025, 0.0125],
                    [0.0125, 0.0125, 0.025],
                ]
            ),
            id="degree 1, first cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[0.3, 0.0], [0.3, 1.0], [0.0, 1.0]]),
            np.array(
                [
                    [0.025, 0.0125, 0.0125],
                    [0.0125, 0.025, 0.0125],
                    [0.0125, 0.0125, 0.025],
                ]
            ),
            id="degree 1, second cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[0.3, 0.0], [1.0, 0.0], [0.3, 1.0]]),
            np.array(
                [
                    [7 / 120, 7 / 240, 7 / 240],
                    [7 / 240, 7 / 120, 7 / 240],
                    [7 / 240, 7 / 240, 7 / 120],
                ]
            ),
            id="degree 1, third cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[1.0, 0.0], [1.0, 1.0], [0.3, 1.0]]),
            np.array(
                [
                    [7 / 120, 7 / 240, 7 / 240],
                    [7 / 240, 7 / 120, 7 / 240],
                    [7 / 240, 7 / 240, 7 / 120],
                ]
            ),
            id="degree 1, fourth cell",
        ),
        pytest.param(
            2,
            1,
            np.array([[0.0, 0.0], [0.3, 0.0], [0.0, 1.0]]),
            np.array(
                [
                    [
                        5.00000000e-03,
                        -8.33333333e-04,
                        -8.33333333e-04,
                        -3.33333333e-03,
                        1.70219741e-17,
                        1.48535698e-17,
                    ],
                    [
                        -8.33333333e-04,
                        5.00000000e-03,
                        -8.33333333e-04,
                        2.51534904e-17,
                        -3.33333333e-03,
                        -2.05998413e-18,
                    ],
                    [
                        -8.33333333e-04,
                        -8.33333333e-04,
                        5.00000000e-03,
                        2.29850861e-17,
                        -3.84891771e-18,
                        -3.33333333e-03,
                    ],
                    [
                        -3.33333333e-03,
                        2.55871713e-17,
                        2.27682456e-17,
                        2.66666667e-02,
                        1.33333333e-02,
                        1.33333333e-02,
                    ],
                    [
                        1.70219741e-17,
                        -3.33333333e-03,
                        -3.74049750e-18,
                        1.33333333e-02,
                        2.66666667e-02,
                        1.33333333e-02,
                    ],
                    [
                        1.48535698e-17,
                        -1.84314369e-18,
                        -3.33333333e-03,
                        1.33333333e-02,
                        1.33333333e-02,
                        2.66666667e-02,
                    ],
                ]
            ),
            id="degree 2, first cell",
        ),
        pytest.param(
            2,
            1,
            np.array([[1.0, 0.0], [1.0, 1.0], [0.3, 1.0]]),
            np.array(
                [
                    [
                        1.16666667e-02,
                        -1.94444444e-03,
                        -1.94444444e-03,
                        -7.77777778e-03,
                        4.05491613e-17,
                        3.53449908e-17,
                    ],
                    [
                        -1.94444444e-03,
                        1.16666667e-02,
                        -1.94444444e-03,
                        6.11490025e-17,
                        -7.77777778e-03,
                        -3.25260652e-18,
                    ],
                    [
                        -1.94444444e-03,
                        -1.94444444e-03,
                        1.16666667e-02,
                        5.29090660e-17,
                        -7.26415456e-18,
                        -7.77777778e-03,
                    ],
                    [
                        -7.77777778e-03,
                        6.02816408e-17,
                        5.29090660e-17,
                        6.22222222e-02,
                        3.11111111e-02,
                        3.11111111e-02,
                    ],
                    [
                        4.09828421e-17,
                        -7.77777778e-03,
                        -7.37257477e-18,
                        3.11111111e-02,
                        6.22222222e-02,
                        3.11111111e-02,
                    ],
                    [
                        3.55618313e-17,
                        -4.11996826e-18,
                        -7.77777778e-03,
                        3.11111111e-02,
                        3.11111111e-02,
                        6.22222222e-02,
                    ],
                ]
            ),
            id="degree 2, fourth cell",
        ),
    ],
)
def test_mass_matrix(
    coords,
    element_degree,
    geometry_degree,
    expected_local_matrix,
    code_dir,
    generate_tabulate_function,
):
    """Test code generation for a mass matrix."""
    e = element("Lagrange", "triangle", element_degree)
    space = function_space(
        coordinate_element(element("Lagrange", "triangle", geometry_degree, shape=(2,))), e
    )
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    tabulate = generate_tabulate_function(form, code_dir)

    mat = np.zeros((e.dim, e.dim))
    tabulate(mat, coords)
    assert np.allclose(mat, expected_local_matrix)


@pytest.mark.parametrize(
    ("element_degree", "geometry_degree", "coords", "expected_local_matrix"),
    [
        pytest.param(
            1,
            1,
            np.array([[0.0, 0.0], [0.3, 0.0], [0.0, 1.0]]),
            np.array(
                [
                    [1.81666666666667, -1.66666666666667, -0.15],
                    [-1.66666666666667, 1.66666666666667, 0],
                    [-0.15, 0, 0.15],
                ]
            ),
            id="degree 1, first cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[0.3, 0.0], [0.3, 1.0], [0.0, 1.0]]),
            np.array(
                [
                    [0.15, -0.15, 0],
                    [-0.15, 1.81666666666667, -1.66666666666667],
                    [0, -1.66666666666667, 1.66666666666667],
                ]
            ),
            id="degree 1, second cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[0.3, 0.0], [1.0, 0.0], [0.3, 1.0]]),
            np.array(
                [
                    [1.06428571428571, -0.714285714285714, -0.35],
                    [-0.714285714285714, 0.714285714285714, 0],
                    [-0.35, 0, 0.35],
                ]
            ),
            id="degree 1, third cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[1.0, 0.0], [1.0, 1.0], [0.3, 1.0]]),
            np.array(
                [
                    [0.35, -0.35, 0],
                    [-0.35, 1.06428571428571, -0.714285714285714],
                    [0, -0.714285714285714, 0.714285714285714],
                ]
            ),
            id="degree 1, fourth cell",
        ),
        pytest.param(
            2,
            1,
            np.array([[0.0, 0.0], [0.3, 0.0], [0.0, 1.0]]),
            np.array(
                [
                    [
                        1.816666666666667,
                        0.5555555555555555,
                        0.04999999999999996,
                        -1.1726730697603216e-15,
                        -0.1999999999999985,
                        -2.2222222222222223,
                    ],
                    [
                        0.5555555555555554,
                        1.6666666666666672,
                        0.0,
                        9.992007221626409e-16,
                        -5.551115123125783e-16,
                        -2.2222222222222214,
                    ],
                    [
                        0.049999999999999975,
                        0.0,
                        0.14999999999999994,
                        -3.469446951953614e-17,
                        -0.19999999999999996,
                        6.938893903907228e-17,
                    ],
                    [
                        -1.1796119636642288e-15,
                        8.881784197001252e-16,
                        -2.7755575615628914e-17,
                        4.844444444444447,
                        -4.444444444444446,
                        -0.39999999999999997,
                    ],
                    [
                        -0.1999999999999985,
                        -5.551115123125783e-16,
                        -0.19999999999999996,
                        -4.444444444444446,
                        4.844444444444447,
                        -9.992007221626409e-16,
                    ],
                    [
                        -2.2222222222222223,
                        -2.2222222222222214,
                        7.632783294297951e-17,
                        -0.4,
                        -9.992007221626409e-16,
                        4.844444444444442,
                    ],
                ]
            ),
            id="degree 2, first cell",
        ),
        pytest.param(
            2,
            1,
            np.array([[1.0, 0.0], [1.0, 1.0], [0.3, 1.0]]),
            np.array(
                [
                    [
                        0.3500000000000001,
                        0.1166666666666667,
                        1.457167719820518e-16,
                        -4.163336342344337e-16,
                        3.400058012914542e-16,
                        -0.4666666666666668,
                    ],
                    [
                        0.1166666666666667,
                        1.0642857142857147,
                        0.23809523809523803,
                        -0.9523809523809524,
                        -2.220446049250313e-16,
                        -0.4666666666666657,
                    ],
                    [
                        1.491862189340054e-16,
                        0.23809523809523803,
                        0.7142857142857141,
                        -0.9523809523809527,
                        9.71445146547012e-17,
                        3.608224830031759e-16,
                    ],
                    [
                        -4.163336342344337e-16,
                        -0.9523809523809526,
                        -0.9523809523809527,
                        2.8380952380952387,
                        -0.9333333333333333,
                        -6.661338147750939e-16,
                    ],
                    [
                        3.191891195797325e-16,
                        -2.7755575615628914e-16,
                        9.71445146547012e-17,
                        -0.9333333333333331,
                        2.8380952380952382,
                        -1.9047619047619058,
                    ],
                    [
                        -0.4666666666666667,
                        -0.4666666666666657,
                        3.3306690738754696e-16,
                        -6.661338147750939e-16,
                        -1.9047619047619058,
                        2.8380952380952387,
                    ],
                ]
            ),
            id="degree 2, fourth cell",
        ),
    ],
)
def test_stiffness_matrix(
    coords,
    element_degree,
    geometry_degree,
    expected_local_matrix,
    code_dir,
    generate_tabulate_function,
):
    """Test code generation for a stiffness matrix."""
    e = element("Lagrange", "triangle", element_degree)
    space = function_space(
        coordinate_element(element("Lagrange", "triangle", geometry_degree, shape=(2,))), e
    )
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(grad(u), grad(v)) * dx

    tabulate = generate_tabulate_function(form, code_dir)

    mat = np.zeros((e.dim, e.dim))
    tabulate(mat, coords)
    assert np.allclose(mat, expected_local_matrix)


@pytest.mark.parametrize(
    ("element_degree", "geometry_degree", "coords", "expected_local_vector"),
    [
        pytest.param(
            1,
            1,
            np.array([[0.0, 0.0], [0.3, 0.0], [0.0, 1.0]]),
            np.array([0.0037500000000000033, 0.0075, 0.0037500000000000033]),
            id="degree 1, first cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[0.3, 0.0], [0.3, 1.0], [0.0, 1.0]]),
            np.array([0.011249999999999982, 0.01125, 0.0075]),
            id="degree 1, second cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[0.3, 0.0], [1.0, 0.0], [0.3, 1.0]]),
            np.array([0.05541666666666667, 0.07583333333333332, 0.05541666666666664]),
            id="degree 1, third cell",
        ),
        pytest.param(
            1,
            1,
            np.array([[1.0, 0.0], [1.0, 1.0], [0.3, 1.0]]),
            np.array([0.09625000000000011, 0.09624999999999997, 0.07583333333333334]),
            id="degree 1, fourth cell",
        ),
        pytest.param(
            2,
            1,
            np.array([[0.0, 0.0], [0.3, 0.0], [0.0, 1.0]]),
            np.array(
                [
                    -0.0007500000000000004,
                    0.001499999999999995,
                    -0.0007499999999999996,
                    0.006000000000000002,
                    0.0029999999999999975,
                    0.006000000000000003,
                ]
            ),
            id="degree 2, first cell",
        ),
        pytest.param(
            2,
            1,
            np.array([[1.0, 0.0], [1.0, 1.0], [0.3, 1.0]]),
            np.array(
                [
                    0.004083333333333292,
                    0.004083333333333303,
                    -0.008166666666666676,
                    0.084,
                    0.084,
                    0.10033333333333339,
                ]
            ),
            id="degree 2, fourth cell",
        ),
    ],
)
def test_linear_form(
    coords,
    element_degree,
    geometry_degree,
    expected_local_vector,
    code_dir,
    generate_tabulate_function,
):
    """Test code generation for a linear form."""
    e = element("Lagrange", "triangle", element_degree)
    space = function_space(
        coordinate_element(element("Lagrange", "triangle", geometry_degree, shape=(2,))), e
    )
    v = TestFunction(space)
    x = SpatialCoordinate(2)
    form = x[0] * v * dx

    tabulate = generate_tabulate_function(form, code_dir)

    vec = np.zeros(e.dim)
    tabulate(vec, coords)
    assert np.allclose(vec, expected_local_vector)
