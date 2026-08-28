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
from uflx.geometry import CoordinateDofComponent, expand_geometry
from uflx.graphs import Graph, GraphNode, NodeOrder, generate_graph
from uflx.integrals import dx
from uflx.maps import apply_push_forwards
from uflx_codegeneration.algorithms import expand_inner_products, tabulate_finite_elements
from uflx_codegeneration.generate import (
    extract_domain,
    integrals_to_quadrature,
    pull_back_to_reference,
    tabulate_quadrature,
)
from uflx_codegeneration.nodes import AddToLocalTensor, ArrayEntry, Loop
from uflx_codegeneration.quadrature import QuadratureLoop, quadrature_rule


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
    if len(domain.elements) != 1:
        raise NotImplementedError("Only domains with exactly one element are supported.")
    (coord_element,) = domain.elements
    return coord_element.dim, domain.geometric_dimension


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
    rules = {dx: quadrature_rule(points, weights)}

    graph = form.graph
    assert graph.is_dag()

    graph = integrals_to_quadrature(graph, rules)
    graph = pull_back_to_reference(graph)
    graph = apply_push_forwards(graph)
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
