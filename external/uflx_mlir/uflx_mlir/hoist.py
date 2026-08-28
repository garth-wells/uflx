"""Shared, MLIR-independent loop-hoisting analysis for emit.py's generator.

Motivation: uflx_codegeneration's own pipeline lowers a form like
inner(grad(u), grad(v))*dx into a nested Loop(i) -> Loop(j) ->
QuadratureLoop(q) -> AddToLocalTensor chain, with the (i,j) loops OUTSIDE
the quadrature loop. Emitting code by walking straight from the
AddToLocalTensor's body expression down to its leaves every time places
everything entirely at the innermost (q) loop level -- including
subexpressions (most importantly: the cell's Jacobian/detJ/cofactor
terms, which depend only on the coordinates argument, not on i, j, or q at
all) that are IDENTICAL across every one of the ndofs^2 * nquadrature_points
iterations. On a P3 tetrahedron stiffness kernel this redundant
recomputation made the naive kernel run ~27x slower per call than FFCx's
own compiled kernel; after the reordering and hoisting this module
computes, the same kernel's measured per-call time (via
`ExecutionEngine.lookup`'s direct calling convention, on real MLIR JIT'd
code) dropped from ~140 us/call to ~2.8 us/call -- a ~50x reduction,
putting it ~1.8x FASTER than FFCx per call.

This module computes, for every node in the AddToLocalTensor body's
expression DAG, the shallowest loop level at which it's legal to compute
it once and reuse the result for every iteration below that level --
manual loop-invariant code motion, computed directly from the graph rather
than relying on the LLVM/MLIR optimizer to rediscover it after the fact
(that was tried first -- running `canonicalize,cse,loop-invariant-code-
motion` ahead of dialect conversion -- and confirmed NOT to work: MLIR's
generic loop-invariant-code-motion pass doesn't hoist memref.load-derived
computation without stronger alias guarantees than are available post
conversion from this generic pipeline).
"""

from __future__ import annotations

from uflx.geometry import CoordinateDofComponent
from uflx.graphs import GraphNode, NodeOrder, generate_graph
from uflx_codegeneration.nodes import AddToLocalTensor, ArrayEntry, Loop
from uflx_codegeneration.quadrature import QuadratureLoop


def walk_loop_chain(root: GraphNode) -> tuple[list[tuple[GraphNode, str]], AddToLocalTensor]:
    """Walk down a chain of nested Loop/QuadratureLoop nodes to the AddToLocalTensor at the bottom.

    Args:
        root: The lowered graph's root node.

    Returns:
        A tuple (chain, add_node), where chain is ordered OUTERMOST to
        INNERMOST as [(loop_node, loop_variable_name), ...].

    Raises:
        NotImplementedError: If the chain doesn't end in a single
            AddToLocalTensor (e.g. a form with more than one integral
            term) -- the same restriction emit.generate_mlir_module
            already documents and enforces itself.
    """
    chain: list[tuple[GraphNode, str]] = []
    node: GraphNode = root
    while isinstance(node, (Loop, QuadratureLoop)):
        chain.append((node, node.variable))
        node = node.body
    if not isinstance(node, AddToLocalTensor):
        raise NotImplementedError(
            "Expected a chain of Loop/QuadratureLoop nodes ending in a single "
            f"AddToLocalTensor, got {type(node)} at the bottom -- this likely means "
            "the form has more than one integral term, which isn't handled yet."
        )
    return chain, node


def _direct_loop_deps(node: GraphNode, loop_vars: set[str]) -> frozenset[str]:
    """Return the loop variables `node` references directly through its own index attributes.

    As opposed to through its graph successors, which compute_levels()
    below folds in separately. Only ArrayEntry and CoordinateDofComponent
    carry index tuples that can name a loop variable (a plain int, or the
    digit-string convention uflx_codegeneration also uses for a size-1
    axis, never do) -- see emit.py's `_OpCtx.resolve_index` for the same
    int-vs-loop-variable-name distinction made at emission time; this must
    classify identically or a node could get hoisted above a loop it
    actually depends on.
    """
    if isinstance(node, ArrayEntry):
        indices = node.index
    elif isinstance(node, CoordinateDofComponent):
        indices = (node._point, node._component)
    else:
        return frozenset()

    return frozenset(
        i for i in indices if isinstance(i, str) and not i.lstrip("-").isdigit() and i in loop_vars
    )


def compute_levels(add_node: AddToLocalTensor, loop_vars: list[str]) -> dict[GraphNode, int]:
    """Compute the shallowest legal loop depth for every node in add_node's body.

    For every node in add_node.body's expression DAG, compute the depth
    (0..len(loop_vars)) at which it's legal to compute it once and reuse
    it below: depth D means "can be computed right after entering
    loop_vars[D-1] but before loop_vars[D]" (depth 0 means "before
    entering any loop at all"; depth len(loop_vars) means "only valid at
    the innermost level", i.e. no hoisting benefit).

    A node's depth is the size of the shortest PREFIX of loop_vars that
    covers every loop variable the node's value actually depends on
    (directly, or transitively through its graph successors) -- loop_vars
    must be in the same outer-to-inner order as the actual nesting
    (walk_loop_chain's `chain`), since a node depending on loop_vars[2]
    alone still can't be computed before loop_vars[0]/[1] are in scope
    textually, even though its VALUE doesn't vary with them.

    Correctness of the resulting hoisting relies on one invariant, true by
    construction here: depends_on(parent) is always a superset of
    depends_on(child) (it's built as a union), so level(parent) >=
    level(child) always -- a node is never scheduled before a dependency
    it needs.

    Args:
        add_node: The AddToLocalTensor whose body expression to analyze.
        loop_vars: The loop variable names enclosing add_node, outermost
            first (see reorder_quadrature_outermost).

    Returns:
        A mapping from every node in add_node.body's expression DAG to its
        computed depth.
    """
    loop_var_set = set(loop_vars)
    depends_on: dict[GraphNode, frozenset[str]] = {}
    graph = generate_graph(add_node.body)
    for node in graph.ordered_nodes(NodeOrder.leaves_first):
        deps = set(_direct_loop_deps(node, loop_var_set))
        for child in node.successors:
            deps |= depends_on[child]
        depends_on[node] = frozenset(deps)

    levels: dict[GraphNode, int] = {}
    for node, deps in depends_on.items():
        bound: set[str] = set()
        level = None
        for d in range(len(loop_vars) + 1):
            if deps <= bound:
                level = d
                break
            if d < len(loop_vars):
                bound.add(loop_vars[d])
        if level is None:
            raise AssertionError(
                f"{node!r} depends on loop variable(s) {deps} not covered by {loop_vars} -- "
                "this means some node references a loop variable that isn't actually one of "
                "the loops enclosing it, which should be impossible."
            )
        levels[node] = level
    return levels


def reorder_quadrature_outermost(
    chain: list[tuple[GraphNode, str]],
) -> list[tuple[GraphNode, str]]:
    """Reorder a loop chain so every QuadratureLoop comes before every (dof) Loop.

    Relative order within each group is preserved.

    Why: uflx_codegeneration's own pipeline happens to nest the dof loops
    (Loop, one per tensor axis -- e.g. test/trial function indices)
    OUTSIDE the quadrature loop. That's the worst order for this kind of
    assembly: almost everything that varies with the quadrature point
    (geometry -- Jacobian/detJ/cofactors -- and any per-dof basis
    value/gradient lookup) ends up needing to be recomputed once per
    (dof, dof, ..., quadrature-point) COMBINATION instead of once per
    quadrature point, because a node can only be hoisted as far out as the
    loops that are already open -- and the quadrature loop being
    innermost means anything touching it is stuck at the very bottom
    regardless of whether it also touches the dof indices. Every real FEM
    assembly kernel loops the quadrature point OUTERMOST for exactly this
    reason.

    Reordering is valid here specifically because every loop in this chain
    is a simple rectangular loop with integer-constant bounds (Loop.start/
    end are asserted to be plain ints elsewhere, e.g.
    lowering.collect_int_constants) and the only side effect anywhere in
    the nest is AddToLocalTensor's additive accumulation into the output
    tensor -- no loop's bounds depend on another loop's variable, and
    there is nothing else in the nest to reorder around. Floating-point
    addition isn't perfectly associative, so this can change the
    accumulated sum's last few bits versus the original iteration order,
    but not by more than ordinary quadrature-order-dependent rounding
    already implies -- emit.generate_mlir_module's own tests validate to
    rtol=1e-9, many orders of magnitude looser than that effect.

    Args:
        chain: A loop chain as returned by walk_loop_chain, outermost first.

    Returns:
        The same loops, reordered with every QuadratureLoop first.
    """
    quad = [(n, v) for n, v in chain if isinstance(n, QuadratureLoop)]
    dof = [(n, v) for n, v in chain if not isinstance(n, QuadratureLoop)]
    return quad + dof


def topo_order(add_node: AddToLocalTensor) -> list[GraphNode]:
    """Return add_node.body's expression DAG nodes, children before parents.

    This is the order emit.py's depth-driven driver (`_build_nest`) walks
    in. A node's children are always emitted earlier in this order than
    the node itself (topological sort), and -- combined with the level
    monotonicity invariant above -- earlier than or at the same depth as
    the node, which is exactly what makes single-pass depth-driven
    emission correct without any recursion.

    Args:
        add_node: The AddToLocalTensor whose body expression to order.

    Returns:
        The nodes of add_node.body's expression DAG, children before parents.
    """
    return list(generate_graph(add_node.body).ordered_nodes(NodeOrder.leaves_first))
