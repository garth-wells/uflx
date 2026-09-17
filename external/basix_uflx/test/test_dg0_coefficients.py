"""Verify real DG0 elements expose cellwise-constant coefficient semantics."""

import numpy as np
import pytest
from uflx import Coefficient, coordinate_element, function_space, grad
from uflx.graphs import as_graph

from basix_uflx import element


@pytest.mark.parametrize("cell,dim", [("triangle", 2), ("tetrahedron", 3)])
@pytest.mark.parametrize("geometry_degree", [1, 2])
@pytest.mark.parametrize("blocked", [False, True])
def test_basix_dg0(cell, dim, geometry_degree, blocked):
    """Basix DG0 has constant values and core removes its coefficient gradient."""
    scalar = element("Lagrange", cell, 0, discontinuous=True)
    points = np.array([[0.1] * dim, [0.2] * dim])
    table = np.asarray(scalar.tabulate(1, points))
    np.testing.assert_allclose(table[0], 1)
    np.testing.assert_allclose(table[1:], 0)
    field = element("Lagrange", cell, 0, shape=(dim,), discontinuous=True) if blocked else scalar
    domain = coordinate_element(element("Lagrange", cell, geometry_degree, shape=(dim,)))
    c = Coefficient(function_space(domain, field))
    assert c.is_cellwise_constant
    result = grad(c)
    assert result.value_shape == ((dim, dim) if blocked else (dim,))
    assert not any(isinstance(n, Coefficient) for n in as_graph(result))
