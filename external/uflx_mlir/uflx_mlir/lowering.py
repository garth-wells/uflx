"""Graph-lowering helpers shared by hoist.py's analysis and emit.py's generator.

This deliberately reuses uflx_codegeneration's existing, generic pipeline
for the actual finite-element math (quadrature selection/tabulation,
geometry expansion, pull-back/push-forward, inner-product expansion) --
the same functions imported below from uflx_codegeneration.generate -- and
only replaces uflx_codegeneration.c's C-string-building final step with an
MLIR generator (see emit.py). One deliberate simplification relative to
the C backend: the quadrature rule is chosen generically via
basix.make_quadrature for the integral's actual cell and a computed
degree, instead of uflx_codegeneration.generate.generate()'s hardcoded
triangle/degree-10 rule (an explicit `# TODO` in their own source, not a
structural limitation of the pipeline).
"""

from __future__ import annotations

import basix
import numpy as np
from uflx.complex import take_real_part
from uflx.domains import AbstractCoordinateElement
from uflx.expressions import AbstractExpression
from uflx.geometry import (
    CoordinateDofComponent,
    Jacobian,
    JacobianDeterminant,
    JacobianInverse,
    JacobianInverseTranspose,
    JacobianTranspose,
    expand_geometry,
)
from uflx.graphs import Graph, GraphNode, NodeOrder, generate_graph
from uflx.graphs.algorithms import replace
from uflx.integrals import AbstractMeasure, dx
from uflx.maps import apply_push_forwards
from uflx.tensors import Matrix
from uflx_codegeneration.algorithms import expand_inner_products, tabulate_finite_elements
from uflx_codegeneration.generate import (
    extract_domain,
    integrals_to_quadrature,
    pull_back_to_reference,
    tabulate_quadrature,
)
from uflx_codegeneration.nodes import AddToLocalTensor, ArrayEntry, Loop
from uflx_codegeneration.quadrature import QuadratureLoop, QuadratureRule, quadrature_rule


def coordinate_shape(form) -> tuple[int, int]:
    """Return the coordinate element's (n_coordinate_dofs, tdim) for a form's domain.

    Read directly off the coordinate element before any lowering, so this
    works regardless of what the lowered graph happens to still mention.

    Args:
        form: A UFLx form (as returned by e.g. `inner(...) * dx`).

    Returns:
        A tuple (n_coordinate_dofs, tdim).
    """
    graph = form.graph
    domain = extract_domain(graph, graph.root)
    if not isinstance(domain, AbstractCoordinateElement):
        raise NotImplementedError("Only domains with a coordinate element are supported.")
    if len(domain.elements) != 1:
        raise NotImplementedError("Only domains with exactly one element are supported.")
    (coord_element,) = domain.elements
    gdim = domain.geometric_dimension
    if coord_element.dim % gdim != 0:
        raise ValueError("Coordinate element dimension is not divisible by its value size.")
    return coord_element.dim // gdim, gdim


def _affine_tetrahedron_jacobian() -> Matrix:
    """Build the Jacobian of a P1 tetrahedron directly from its vertices."""
    return Matrix(
        [
            [
                CoordinateDofComponent(column + 1, row, 3) - CoordinateDofComponent(0, row, 3)
                for column in range(3)
            ]
            for row in range(3)
        ]
    )


def _expand_affine_tetrahedron_geometry(graph: Graph, cell: basix.CellType) -> Graph:
    """Expand P1 tetrahedral Jacobian nodes without a quadrature-dependent FE table."""
    if cell != basix.CellType.tetrahedron:
        return graph

    replacements: dict[GraphNode, GraphNode] = {}
    for node in graph:
        if not isinstance(
            node,
            (
                Jacobian,
                JacobianDeterminant,
                JacobianInverse,
                JacobianTranspose,
                JacobianInverseTranspose,
            ),
        ):
            continue
        if len(node.domain.elements) != 1:
            continue
        (coordinate_element,) = node.domain.elements
        if coordinate_element.lagrange_superdegree != 1 or node.domain.geometric_dimension != 3:
            continue

        jacobian = _affine_tetrahedron_jacobian()
        replacement: AbstractExpression
        if isinstance(node, Jacobian):
            replacement = jacobian
        elif isinstance(node, JacobianDeterminant):
            replacement = abs(jacobian.compute_determinant())
        elif isinstance(node, JacobianInverse):
            replacement = jacobian.compute_inverse()
        elif isinstance(node, JacobianTranspose):
            replacement = jacobian.transpose()
        else:
            assert isinstance(node, JacobianInverseTranspose)
            replacement = jacobian.compute_inverse().transpose()
        replacements[node] = replacement

    return replace(graph, replacements) if replacements else graph


def lower_form(form, degree: int, cell: basix.CellType) -> tuple[dict[str, np.ndarray], Graph]:
    """Lower a UFLx form to a graph of arithmetic/loop/table nodes.

    Runs uflx_codegeneration's own pipeline stages, with a quadrature rule
    chosen for the actual cell instead of a hardcoded one.

    Args:
        form: A UFLx form (as returned by e.g. `inner(...) * dx`).
        degree: The polynomial degree used to size the quadrature rule.
        cell: The reference cell the form is integrated over.

    Returns:
        A tuple (tables, graph): `tables` maps table name to its constant
        numpy array (quadrature weights, tabulated basis functions, ...);
        `graph` is the fully-lowered DAG whose root is a chain of
        Loop/QuadratureLoop nodes ending in a single AddToLocalTensor.
    """
    qdeg = max(2 * (degree - 1), 1)
    points, weights = basix.make_quadrature(cell, qdeg)
    rules: dict[AbstractMeasure, QuadratureRule] = {dx: quadrature_rule(points, weights)}

    graph = form.graph
    assert graph.is_dag()

    # pull_back_to_reference/apply_push_forwards must run BEFORE
    # integrals_to_quadrature: since uflx#51 ("Move mapping of integral to
    # reference into UFLx core"), Jacobian/JacobianDeterminant/
    # JacobianInverse/JacobianTranspose/JacobianInverseTranspose are built
    # with point=None by pull_back_to_reference, and it's
    # integrals_to_quadrature that walks an integral's descendants and
    # replaces each point-less Jacobian-family node with one bound to the
    # actual quadrature point -- matching the order
    # uflx_codegeneration.generate.generate() itself now uses. Running
    # integrals_to_quadrature first (the old order) leaves those nodes'
    # `point` permanently None, and expand_geometry's
    # `assert self.point is not None` fires downstream.
    graph = pull_back_to_reference(graph)
    graph = apply_push_forwards(graph)
    graph = integrals_to_quadrature(graph, rules)
    graph = _expand_affine_tetrahedron_geometry(graph, cell)
    graph = expand_geometry(graph)
    graph = expand_inner_products(graph)
    graph = take_real_part(graph)

    q_tables, graph = tabulate_quadrature(graph)
    fe_tables, graph = tabulate_finite_elements(graph)
    tables = {**q_tables, **fe_tables}
    return tables, graph


def collect_int_constants(root: GraphNode, a_shape: tuple[int, ...]) -> set[int]:
    """Collect every literal int the lowered graph needs as an `index` constant.

    Covers loop bounds and int entries in ArrayEntry/AddToLocalTensor/
    CoordinateDofComponent index tuples.

    Args:
        root: The lowered graph's root node (a Loop/QuadratureLoop chain).
        a_shape: The shape of the local tensor being assembled.

    Returns:
        The set of distinct non-negative ints that appear as index/loop-bound
        literals anywhere in the graph, always including 0, 1, and every
        entry of `a_shape`.
    """
    ints: set[int] = {0, 1, *a_shape}
    graph = generate_graph(root)
    for node in graph.ordered_nodes(NodeOrder.roots_first):
        if isinstance(node, Loop):
            assert isinstance(node.start, int) and isinstance(node.end, int)
            ints.add(node.start)
            ints.add(node.end)
        elif isinstance(node, QuadratureLoop):
            ints.add(node.rule.npoints)
        elif isinstance(node, AddToLocalTensor):
            for i in node.component:
                if isinstance(i, int):
                    ints.add(i)
                elif isinstance(i, str) and i.lstrip("-").isdigit():
                    ints.add(int(i))
        elif isinstance(node, ArrayEntry):
            for i in node.index:
                if isinstance(i, int):
                    ints.add(i)
                elif isinstance(i, str) and i.lstrip("-").isdigit():
                    ints.add(int(i))
        elif isinstance(node, CoordinateDofComponent):
            ints.add(node._point)
            ints.add(node._component)
    return ints
