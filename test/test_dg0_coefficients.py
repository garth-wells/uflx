"""Cellwise-constant coefficient semantics, independently of code generation."""

from itertools import product

import pytest
from conftest import LagrangeElement

from uflx import Coefficient, TestFunction, coordinate_element, dx, function_space, grad, inner
from uflx.algorithms import pull_back_to_reference, reconstruct_node, replace
from uflx.expressions import AbstractExpression
from uflx.graphs import as_graph
from uflx.integrals import Integral
from uflx.maps import AbstractReferenceMap, BlockedReferenceMap, IdentityReferenceMap
from uflx.operators import Grad, ReferenceGrad
from uflx.tensors import zero


def assert_zero(expression, shape):
    """Check every scalar component and ensure no coefficient remains."""
    assert expression.value_shape == shape
    for index in product(*(range(n) for n in shape)):
        value = expression.component(*index) if shape else expression
        assert value.as_float() == 0
    assert not any(isinstance(node, Coefficient) for node in as_graph(expression))


@pytest.mark.parametrize(
    "cell,dim", [("triangle", 2), ("tetrahedron", 3), ("quadrilateral", 2), ("hexahedron", 3)]
)
@pytest.mark.parametrize("geometry_degree", [1, 2])
def test_dg0_value_and_gradient(lagrange_element, cell, dim, geometry_degree):
    """DG0 values survive in forms, but their gradients contain no coefficient."""
    domain = coordinate_element(lagrange_element(cell, geometry_degree, (dim,)))
    space = function_space(domain, lagrange_element(cell, 0))
    c = Coefficient(space)
    other = Coefficient(space)
    v = TestFunction(function_space(domain, lagrange_element(cell, 1)))
    assert c.is_cellwise_constant
    assert c.label != other.label
    form = inner(c, v) * dx
    assert isinstance(form, Integral)
    for expression in [form, pull_back_to_reference(form)]:
        coefficients = [n for n in as_graph(expression) if isinstance(n, Coefficient)]
        assert len(coefficients) == 1
        assert coefficients[0].label == c.label
        assert coefficients[0].is_cellwise_constant
        assert coefficients[0].integral_label == form.label
    reconstructed = reconstruct_node(c, {})
    assert isinstance(reconstructed, Coefficient)
    assert reconstructed.is_cellwise_constant
    assert_zero(grad(c), (dim,))
    assert_zero(pull_back_to_reference(Grad(c)), (dim,))
    reference = Coefficient(space, is_reference=True)
    assert_zero(ReferenceGrad(reference).expand_geometry(), (dim,))
    for i in range(dim):
        assert_zero(c.diff(i), ())
        assert_zero(reference.diff(i), ())
    gradient_form = inner(grad(c), grad(v)) * dx
    assert not any(isinstance(n, Coefficient) for n in as_graph(gradient_form))
    # Replacement must retain the new coefficient's constancy on reconstruction.
    variable = Coefficient(v.function_space)
    replaced = replace(Grad(variable), {variable: c})
    assert_zero(pull_back_to_reference(replaced), (dim,))
    with pytest.raises(ValueError):
        c.diff(dim)


@pytest.mark.parametrize("degree", [1, 2])
def test_nonconstant_control(lagrange_element, degree):
    """Higher-order coefficients must retain their gradients."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    c = Coefficient(function_space(domain, lagrange_element("triangle", degree)))
    assert not c.is_cellwise_constant
    assert isinstance(grad(c), Grad)
    with pytest.raises(NotImplementedError):
        c.diff(0)


class UnknownMap(AbstractReferenceMap):
    """A mapping with no guarantee that it preserves reference constants."""

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Leave this mock map unimplemented."""
        raise NotImplementedError()

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Leave this mock map unimplemented."""
        raise NotImplementedError()

    def physical_value_shape(self, geometric_dimension: int) -> tuple[int, ...]:
        """Use scalar values."""
        return ()


def test_unknown_mapping_is_not_constant(lagrange_element):
    """Reference degree zero alone is insufficient on the physical cell."""
    base = lagrange_element("triangle", 0)

    class MappedElement(LagrangeElement):
        @property
        def reference_map(self):
            return UnknownMap()

    domain = coordinate_element(lagrange_element("triangle", 2, (2,)))
    space = function_space(domain, MappedElement(base.cell, 0))
    c = Coefficient(space)
    assert not c.is_cellwise_constant
    assert isinstance(grad(c), Grad)
    reference = Coefficient(space, is_reference=True)
    assert reference.is_cellwise_constant
    assert_zero(ReferenceGrad(reference).expand_geometry(), (2,))
    assert not BlockedReferenceMap(UnknownMap(), (2,)).preserves_constant_values
    assert BlockedReferenceMap(IdentityReferenceMap(), (2,)).preserves_constant_values


@pytest.mark.parametrize("shape", [(), (2,), (2, 3), (2, 2, 2)])
def test_shaped_zero(shape):
    """Zero gradients and derivatives preserve all value axes."""
    assert_zero(zero(shape), shape)


@pytest.mark.parametrize("shape", [(2,), (2, 2)])
def test_blocked_constant_derivatives(lagrange_element, shape):
    """Vector and tensor constants produce zeros with the appropriate ranks."""

    class BlockedElement(LagrangeElement):
        @property
        def reference_map(self):
            return BlockedReferenceMap(IdentityReferenceMap(), shape)

    scalar = lagrange_element("triangle", 0)
    domain = coordinate_element(lagrange_element("triangle", 2, (2,)))
    space = function_space(domain, BlockedElement(scalar.cell, 0, shape))
    c = Coefficient(space)
    assert c.is_cellwise_constant
    assert_zero(grad(c), (*shape, 2))
    assert_zero(pull_back_to_reference(Grad(c)), (*shape, 2))
    assert_zero(c.diff(0), shape)
    reference = Coefficient(space, is_reference=True)
    assert_zero(ReferenceGrad(reference).expand_geometry(), (*shape, 2))


def test_unknown_degree_and_multielement_space(lagrange_element):
    """Neither unknown degrees nor a mixture containing P1 are known constant."""

    class UnknownDegreeElement(LagrangeElement):
        @property
        def lagrange_superdegree(self):
            return None

    p0 = lagrange_element("triangle", 0)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    for elements in [UnknownDegreeElement(p0.cell, 0), (p0, lagrange_element("triangle", 1))]:
        c = Coefficient(function_space(domain, elements))
        assert not c.is_cellwise_constant
        assert isinstance(grad(c), Grad)


@pytest.mark.parametrize("shape", [(0,), (-1,), (2, 0)])
def test_zero_rejects_invalid_shape(shape):
    """Reject tensor extents the core tensor representation cannot express."""
    with pytest.raises(ValueError):
        zero(shape)
