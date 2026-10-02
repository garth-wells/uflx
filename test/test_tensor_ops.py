"""Test uflx.operators's sym/skew/tr/transpose/dev tensor operators."""

import numpy as np
import pytest

from uflx import coordinate_element, function_space
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.expressions import Product, RealScalar
from uflx.operators import ReferenceGrad, Tr, dev, skew, sym, tr, transpose
from uflx.points import Point
from uflx.tensors import Identity, Matrix, Vector, zero


def _to_matrix(values: np.ndarray) -> Matrix:
    return Matrix([[RealScalar(float(v)) for v in row] for row in values])


def _to_numpy(expression, shape: tuple[int, ...]) -> np.ndarray:
    if shape == ():
        return np.array(expression.as_float())
    rows, cols = shape
    return np.array(
        [[expression.component(i, j).as_float() for j in range(cols)] for i in range(rows)]
    )


@pytest.mark.parametrize("n", [1, 2, 3])
def test_transpose_matches_numpy(n):
    """transpose() should match numpy's transpose, including for non-square shapes."""
    rng = np.random.default_rng(2)
    values = rng.uniform(-2.0, 2.0, (n, n + 1))
    matrix = _to_matrix(values)

    result = transpose(matrix)
    assert result.value_shape == (n + 1, n)
    np.testing.assert_allclose(_to_numpy(result, result.value_shape), values.T, rtol=1e-12)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_tr_matches_numpy(n):
    """tr() should match numpy's trace for square matrices."""
    rng = np.random.default_rng(3)
    values = rng.uniform(-2.0, 2.0, (n, n))
    matrix = _to_matrix(values)

    np.testing.assert_allclose(tr(matrix).as_float(), np.trace(values), rtol=1e-12)


@pytest.mark.parametrize("n", [2, 3])
def test_sym_and_skew_match_numpy(n):
    """sym()/skew() should match (A + A.T)/2 and (A - A.T)/2, and recombine to A."""
    rng = np.random.default_rng(4)
    values = rng.uniform(-2.0, 2.0, (n, n))
    matrix = _to_matrix(values)

    s = sym(matrix)
    k = skew(matrix)
    assert s.value_shape == k.value_shape == (n, n)
    np.testing.assert_allclose(_to_numpy(s, (n, n)), (values + values.T) / 2, rtol=1e-12)
    np.testing.assert_allclose(_to_numpy(k, (n, n)), (values - values.T) / 2, rtol=1e-12)

    recombined = _to_numpy(s, (n, n)) + _to_numpy(k, (n, n))
    np.testing.assert_allclose(recombined, values, rtol=1e-12)


@pytest.mark.parametrize("n", [2, 3])
def test_dev_matches_numpy(n):
    """dev() should match A - tr(A)/n * I, and its trace should always vanish."""
    rng = np.random.default_rng(5)
    values = rng.uniform(-2.0, 2.0, (n, n))
    matrix = _to_matrix(values)

    d = dev(matrix)
    assert d.value_shape == (n, n)
    expected = values - np.trace(values) / n * np.eye(n)
    np.testing.assert_allclose(_to_numpy(d, (n, n)), expected, rtol=1e-12)
    np.testing.assert_allclose(tr(d).as_float(), 0.0, atol=1e-12)


def test_identity_is_identity_matrix():
    """Identity(d) should be the usual d x d identity tensor."""
    for d in (1, 2, 3):
        identity = Identity(d)
        assert identity.value_shape == (d, d)
        np.testing.assert_allclose(_to_numpy(identity, (d, d)), np.eye(d))


def test_identity_repr():
    """Identity's repr should reflect its own (renamed) class name."""
    assert repr(Identity(3)) == "Identity(3)"


def test_identity_equality():
    """Identity instances should compare equal by size, and only by size."""
    assert Identity(2) == Identity(2)
    assert Identity(2) != Identity(3)


def test_identity_matrix_product_acts_as_identity():
    """Identity(d) should act as a true identity under matrix multiplication."""
    values = np.array([[1.0, 2.0], [3.0, 4.0]])
    matrix = _to_matrix(values)

    np.testing.assert_allclose(_to_numpy(Identity(2) @ matrix, (2, 2)), values)
    np.testing.assert_allclose(_to_numpy(matrix @ Identity(2), (2, 2)), values)


def test_transpose_requires_rank_2():
    """transpose() should reject vectors and scalars, not just silently misbehave."""
    with pytest.raises(ValueError):
        transpose(Vector([RealScalar(1.0), RealScalar(2.0)]))
    with pytest.raises(ValueError):
        transpose(RealScalar(1.0))


def test_tr_and_dev_require_square():
    """tr()/dev() should reject non-square rank-2 expressions."""
    non_square = Matrix([[RealScalar(1.0), RealScalar(2.0)]])
    with pytest.raises(ValueError):
        tr(non_square)
    with pytest.raises(ValueError):
        dev(non_square)


def test_sym_skew_reject_non_square():
    """sym()/skew() of a non-square matrix should fail (shape mismatch in A + A.T)."""
    non_square = Matrix([[RealScalar(1.0), RealScalar(2.0)]])
    with pytest.raises(AssertionError):
        sym(non_square)
    with pytest.raises(AssertionError):
        skew(non_square)


def test_tr_component_raises():
    """A scalar Tr node should refuse .component(), like every other scalar expression."""
    matrix = _to_matrix(np.eye(2))
    with pytest.raises(ValueError):
        Tr(matrix).component(0)


def test_scalar_product_times_matrix_regression():
    """A scalar built from `*` (a Product instance) times a matrix.

    It must still dispatch to ScalarMult rather than being flattened into
    Product's own (same-shape-only) items list.
    """
    scalar = RealScalar(2.0) * RealScalar(3.0)
    assert isinstance(scalar, Product)
    matrix = _to_matrix(np.eye(2) * 5.0)
    result = scalar * matrix
    assert result.value_shape == (2, 2)
    np.testing.assert_allclose(_to_numpy(result, (2, 2)), np.eye(2) * 5.0 * 6.0, rtol=1e-12)


def test_sigma_isotropic_elasticity_matches_numpy():
    """sigma(u) = lambda*tr(sym(grad(u)))*I + 2*mu*sym(grad(u)), on a literal gradient."""
    rng = np.random.default_rng(6)
    grad_u = rng.uniform(-1.0, 1.0, (3, 3))
    grad_u_expr = _to_matrix(grad_u)
    lambda_, mu = 1.7, 0.8

    eps = sym(grad_u_expr)
    sigma = lambda_ * tr(eps) * Identity(3) + 2 * mu * eps

    eps_np = (grad_u + grad_u.T) / 2
    sigma_np = lambda_ * np.trace(eps_np) * np.eye(3) + 2 * mu * eps_np
    np.testing.assert_allclose(_to_numpy(sigma, (3, 3)), sigma_np, rtol=1e-10)


def test_reference_grad_expand_geometry_vector_shape_regression(lagrange_element):
    """ReferenceGrad.expand_geometry() must match shapes for a vector-valued argument.

    Previously it silently claimed (*value_shape, domain_size) but actually
    built a plain (domain_size,) Vector of non-scalar entries.
    """
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, lagrange_element("triangle", 1, (2,)))
    point = Point([RealScalar(0.25), RealScalar(0.25)])
    basis_function = EvaluatedBasisFunction(space, 0, point, True)

    reference_grad = ReferenceGrad(basis_function)
    expanded = reference_grad.expand_geometry()
    assert expanded.value_shape == reference_grad.value_shape == (2, 2)

    # Every entry must itself be a genuine scalar, not a smuggled-in vector.
    for i in range(2):
        for j in range(2):
            assert expanded.component(i, j).value_shape == ()


def test_reference_grad_expand_geometry_scalar_unchanged(lagrange_element):
    """The pre-existing scalar-argument path is untouched by the vector-shape fix."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, lagrange_element("triangle", 1))
    point = Point([RealScalar(0.25), RealScalar(0.25)])
    basis_function = EvaluatedBasisFunction(space, 0, point, True)

    reference_grad = ReferenceGrad(basis_function)
    expanded = reference_grad.expand_geometry()
    assert expanded.value_shape == reference_grad.value_shape == (2,)


def test_reference_grad_cellwise_constant_vector(lagrange_element):
    """A cellwise-constant vector argument's reference gradient is still exactly zero."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, lagrange_element("triangle", 0, (2,)))
    point = Point([RealScalar(0.25), RealScalar(0.25)])
    basis_function = EvaluatedBasisFunction(space, 0, point, True)
    assert basis_function.is_cellwise_constant

    reference_grad = ReferenceGrad(basis_function)
    expanded = reference_grad.expand_geometry()
    assert expanded.value_shape == (2, 2)
    assert expanded == zero((2, 2))
