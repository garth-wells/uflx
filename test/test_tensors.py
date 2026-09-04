"""Test uflx.tensors.Matrix's compute_inverse/compute_determinant."""

import numpy as np
import pytest

from uflx.expressions import Add, Div, Integer, Mult, Neg, RealScalar, Subtract
from uflx.tensors import Matrix


def _evaluate(expr) -> float:
    """Numerically evaluate an expression tree.

    Handles trees built only from RealScalar literals and +, -, *, / --
    enough for these tests, not a general evaluator.
    """
    if isinstance(expr, (RealScalar, Integer)):
        return float(expr.value)
    if isinstance(expr, Add):
        return _evaluate(expr.first) + _evaluate(expr.second)
    if isinstance(expr, Subtract):
        return _evaluate(expr.first) - _evaluate(expr.second)
    if isinstance(expr, Mult):
        return _evaluate(expr.first) * _evaluate(expr.second)
    if isinstance(expr, Div):
        return _evaluate(expr.first) / _evaluate(expr.second)
    if isinstance(expr, Neg):
        return -_evaluate(expr.argument)
    raise NotImplementedError(f"Cannot evaluate {type(expr)}")


def _to_matrix(values: np.ndarray) -> Matrix:
    return Matrix([[RealScalar(float(v)) for v in row] for row in values])


def _to_numpy(matrix: Matrix) -> np.ndarray:
    rows, cols = matrix.value_shape
    return np.array([[_evaluate(matrix.component(i, j)) for j in range(cols)] for i in range(rows)])


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

    det = _evaluate(_to_matrix(values).compute_determinant())

    np.testing.assert_allclose(det, np.linalg.det(values), rtol=1e-12)


def test_compute_inverse_not_implemented_for_4x4():
    """compute_inverse() should still raise a clear error for unsupported sizes."""
    matrix = _to_matrix(np.eye(4))
    with pytest.raises(NotImplementedError):
        matrix.compute_inverse()
