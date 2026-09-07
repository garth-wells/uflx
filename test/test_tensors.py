"""Test uflx.tensors.Matrix's compute_inverse/compute_determinant."""

import numpy as np
import pytest

from uflx.expressions import Add, Div, Integer, Mult, Neg, RealScalar, Subtract
from uflx.tensors import Matrix


def _to_matrix(values: np.ndarray) -> Matrix:
    return Matrix([[RealScalar(float(v)) for v in row] for row in values])


def _to_numpy(matrix: Matrix) -> np.ndarray:
    rows, cols = matrix.value_shape
    return np.array([[matrix.component(i, j).as_float() for j in range(cols)] for i in range(rows)])


@pytest.mark.parametrize("n", [1, 2, 3])
def test_compute_inverse(n):
    """compute_inverse() should match numpy's inverse for 1x1, 2x2 and 3x3."""
    rng = np.random.default_rng(0)
    values = rng.uniform(0.5, 2.0, (n, n))

    inverse = _to_numpy(_to_matrix(values).compute_inverse())

    np.testing.assert_allclose(inverse, np.linalg.inv(values), rtol=1e-12)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_compute_determinant(n):
    """compute_determinant() should match numpy's determinant for 1x1, 2x2 and 3x3."""
    rng = np.random.default_rng(1)
    values = rng.uniform(0.5, 2.0, (n, n))

    det = _to_matrix(values).compute_determinant().as_float()

    np.testing.assert_allclose(det, np.linalg.det(values), rtol=1e-12)


def test_compute_inverse_not_implemented_for_4x4():
    """compute_inverse() should still raise a clear error for unsupported sizes."""
    matrix = _to_matrix(np.eye(4))
    with pytest.raises(NotImplementedError):
        matrix.compute_inverse()
