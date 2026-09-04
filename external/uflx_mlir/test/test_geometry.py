"""Tests for geometry extraction that do not require LLVM's MLIR bindings."""

from __future__ import annotations

import basix
from basix_uflx import element
from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner
from uflx.geometry import CoordinateDofComponent

from uflx_mlir.geometry import GeometryTensorComponent, geometry_kernel_name
from uflx_mlir.hoist import compute_fission_plan, reorder_quadrature_outermost, walk_loop_chain
from uflx_mlir.lowering import lower_form


def _build_form(degree: int, *, stiffness: bool = True):
    e = element("Lagrange", "tetrahedron", degree, lagrange_variant="equispaced")
    domain = coordinate_element(element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = function_space(domain, e)
    u = TrialFunction(space)
    v = TestFunction(space)
    integrand = inner(grad(u), grad(v)) if stiffness else inner(u, v)
    return integrand * dx


def test_affine_poisson_geometry_is_extracted_from_tabulation_graph() -> None:
    """The element graph should consume packed geometry, not coordinate dofs."""
    _, graph, geometry = lower_form(_build_form(2), 2, basix.CellType.tetrahedron)

    assert geometry is not None
    assert geometry.output_size == 6
    assert geometry.scope == "cell"
    components = [node for node in graph if isinstance(node, GeometryTensorComponent)]
    assert {node.packed_index for node in components} == set(range(6))
    assert not any(isinstance(node, CoordinateDofComponent) for node in graph)


def test_geometry_contraction_uses_three_scratch_vectors() -> None:
    """Keep G @ grad factorised so fission scales with vector dimension, not tensor size."""
    _, graph, _ = lower_form(_build_form(2), 2, basix.CellType.tetrahedron)
    chain, add_node = walk_loop_chain(graph.root)
    chain = reorder_quadrature_outermost(chain)
    _, groups = compute_fission_plan(add_node, [variable for _, variable in chain])

    assert sum(len(group.scratch) for group in groups) == 3


def test_unsupported_form_keeps_inline_geometry() -> None:
    """Extraction should leave non-Poisson forms on the existing lowering path."""
    _, graph, geometry = lower_form(_build_form(1, stiffness=False), 1, basix.CellType.tetrahedron)

    assert geometry is None
    assert any(isinstance(node, CoordinateDofComponent) for node in graph)


def test_geometry_kernel_name() -> None:
    """The companion geometry symbol should be predictable for callers."""
    assert geometry_kernel_name("tabulate_tensor") == "tabulate_tensor_geometry"
