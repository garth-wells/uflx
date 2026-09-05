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
putting it ~1.8x FASTER than FFCx per call at the time this was written.
(A later ffcx release -- 0.12.0.dev0 -- closed most of that gap and pulled
slightly ahead again at P3, ~1.107 us/call vs this module's ~1.978 us/call;
diffing its generated C against this module's MLIR output showed FFCx
folding the quadrature weight into its 6-entry geometric factor ONCE PER
QUADRATURE POINT, rather than once per (dof, dof) pair the way this module
was doing -- see distribute_shallow_factors below, added to close exactly
that gap by performing the equivalent algebraic push-down here.)

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

from dataclasses import dataclass
from typing import cast

from uflx.expressions import Add, Mult, Neg, Subtract
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


def distribute_shallow_factors(add_node: AddToLocalTensor, loop_vars: list[str]) -> None:
    """Push a shallow multiplicative factor down into a deeper Mul/Add tree.

    Motivation: uflx_codegeneration's own form lowering places the
    integration measure's quadrature weight as the OUTERMOST multiplicative
    factor over the whole per-entry integrand -- add_node.body is literally
    Mult(weight, big_expression), where `weight` depends only on the
    quadrature loop variable but `big_expression` depends on every loop
    variable (the quadrature point AND both dof indices). compute_levels()
    can only ever assign nodes that already exist in the graph a legal
    depth -- it can't change what's being multiplied by what -- so
    Mult(weight, big_expression)'s own level is the union of both operands'
    dependencies (everything), and it never gets hoisted even though
    `weight` itself could, in principle, be folded into big_expression's
    own quadrature-point-only subexpressions (e.g. a geometric factor)
    instead of being multiplied into the full, already dof-dependent
    result once per (dof, dof) pair.

    This is exactly the algorithmic difference found (by diffing generated
    code) between this module's output and ffcx 0.12.0.dev0's own
    tabulate_tensor codegen at P3: FFCx pre-multiplies its 6-entry
    geometric factor by the quadrature weight ONCE per quadrature point (6
    scalar multiplies total per point), where this module's un-rewritten
    output multiplied the weight in once per (dof, dof) pair (400 multiplies
    per point) -- see hoist.py's module docstring. This function performs
    the equivalent algebraic rewrite: wherever a Mult node has one operand
    at a shallower computed level than the other, the shallower operand
    gets pushed down through the deeper operand's own Add/Subtract/Neg/Mult
    structure (distributing over Add/Subtract, commuting through Neg, and
    at each Mult recursing into whichever of its two operands has the
    smaller level) until it reaches a leaf or a node type it doesn't know
    how to push through (e.g. Div, Abs), where it's finally multiplied in
    directly.

    This is always mathematically exact -- distributing a scalar factor
    over addition/subtraction/negation, and re-associating a chain of
    multiplications, are both exactly the identities they appear to be, not
    approximations -- and can only ever keep or reduce the level
    compute_levels()/compute_fission_plan() will assign to the rewritten
    nodes relative to the original, un-rewritten Mult, since every new node
    built here combines two operands whose levels were already <= that
    original Mult's level. It mutates add_node.body in place; call this
    BEFORE compute_fission_plan(), so the fission/depth analysis that
    follows sees the already-rewritten (and now more hoistable) tree.

    Only ever inspects/rewrites the TOP-level node of add_node.body -- if
    that top node isn't itself a Mult with two differently-leveled
    operands (e.g. a mass-matrix-style u*v*dx form has the same shape;
    something else may not), this is a no-op.

    Floating point note: this changes the ORDER subexpressions are summed
    and multiplied in (e.g. (w*a)+(w*b) instead of w*(a+b)), which can
    shift the last bit or two of the result versus the un-rewritten
    computation -- not a change in which mathematical quantity gets
    computed, just the rounding path to it. The existing rtol=1e-9-or-
    looser reference comparisons this package's tests already run (see
    reorder_quadrature_outermost's own docstring for the same caveat about
    its reordering) comfortably cover a difference this small.

    Args:
        add_node: The AddToLocalTensor whose body to rewrite in place.
        loop_vars: The loop variables enclosing add_node, outermost first
            (see reorder_quadrature_outermost) -- same as compute_levels'
            own loop_vars argument.
    """
    levels = compute_levels(add_node, loop_vars)

    def push(factor, node):
        if isinstance(node, Add):
            return Add(push(factor, node.first), push(factor, node.second))
        if isinstance(node, Subtract):
            return Subtract(push(factor, node.first), push(factor, node.second))
        if isinstance(node, Neg):
            return Neg(push(factor, node.argument))
        if isinstance(node, Mult):
            if levels[node.first] <= levels[node.second]:
                return Mult(push(factor, node.first), node.second)
            return Mult(node.first, push(factor, node.second))
        # Leaf (ArrayEntry, GeometryTensorComponent, CoordinateDofComponent,
        # RealScalar, Integer, ...) or an operator not handled above (Div,
        # Abs, ...) -- multiply here directly. Always correct (see
        # docstring), just not always maximally hoisted.
        return Mult(factor, node)

    body = add_node.body
    if not isinstance(body, Mult):
        return
    first_level = levels[body.first]
    second_level = levels[body.second]
    if first_level < second_level:
        add_node.body = push(body.first, body.second)
    elif second_level < first_level:
        add_node.body = push(body.second, body.first)
    # else: equal levels -- nothing to gain by pushing either way, leave
    # add_node.body exactly as it was.


def _prefix_depth(deps: frozenset[str], loop_vars: list[str]) -> int:
    """The largest d such that set(loop_vars[:d]) is a SUBSET of deps.

    This is the opposite comparison from compute_levels()'s notion of
    level (the smallest d such that deps is a subset of
    set(loop_vars[:d])). The two coincide whenever deps IS some prefix of
    loop_vars exactly (the common case): both then equal that prefix's
    length. They diverge exactly when deps references a loop variable
    "past a gap" -- e.g. deps={q, j} with loop_vars=[q, i, j] (the node
    needs q and j but not i, which sits between them in the nesting):
    compute_levels' level is 3 (deps isn't covered by any prefix shorter
    than the whole chain), while this function returns 1 (every loop
    variable *before* the gap -- here just {q} -- is still a valid landing
    point; d=1 is as far out as this node's computation can be hoisted
    without first resolving the gap). See compute_fission_plan, which uses
    this to detect exactly that situation and route it through an auxiliary
    loop instead of forcing the node down to the innermost level.
    """
    d = 0
    for v in loop_vars:
        if v not in deps:
            break
        d += 1
    return d


@dataclass(frozen=True)
class FissionGroup:
    """A set of expression-DAG nodes sharing the same "gap" dependence pattern.

    E.g. every node here depends on the
    quadrature point and the trial-function dof index, but not the
    test-function dof index sitting between them in the loop nest -- and so
    can't be hoisted any further by compute_levels()'s ordinary
    prefix-based depth alone (see _prefix_depth's docstring). Instead they
    get computed once each, ahead of time, in a small auxiliary loop nest
    covering exactly the "gap" variables, with results saved to a scratch
    buffer indexed by those variables and read back via a plain array
    lookup wherever they're actually needed in the main nest -- the same
    technique the hand-written kernels in mlir-kernels/python/
    generate_kernel.py already use (their numx_scratch/numy_scratch/
    numz_scratch buffers).

    Attributes:
        depth: the main-chain depth at which this group's auxiliary loop
            should be inserted -- i.e. right after loop_vars[depth - 1] is
            entered and before loop_vars[depth] is entered. Every node in
            this group depends on every loop variable in loop_vars[:depth]
            (already in scope there) plus every variable in gap_vars, and
            nothing else.
        gap_vars: the loop variables (in the same outer-to-inner order as
            the main chain) this group depends on beyond `depth` -- the
            auxiliary loop nests one scf.for per entry, in this order.
        nodes: every node assigned to this group, in dependency
            (children-before-parents) order -- exactly the subset of
            topo_order()'s result that must be computed inside the
            auxiliary loop rather than the main nest.
        scratch: the subset of `nodes` that must actually be written to a
            scratch buffer -- i.e. has at least one dependent (successor's
            parent) outside this exact group, so its value must survive
            past the auxiliary loop closing. A node in `nodes` but not in
            `scratch` is a pure local intermediate: every consumer of it is
            also in this same group, so it only ever needs to exist as an
            ordinary SSA value inside the auxiliary loop's own body.
    """

    depth: int
    gap_vars: tuple[str, ...]
    nodes: tuple[GraphNode, ...]
    scratch: tuple[GraphNode, ...]


def compute_fission_plan(
    add_node: AddToLocalTensor, loop_vars: list[str]
) -> tuple[dict[GraphNode, int], list[FissionGroup]]:
    """Extend compute_levels() with loop-fission groups.

    Covers nodes a plain prefix-based depth can't hoist at all (see
    _prefix_depth and FissionGroup above).

    Returns (levels, fission_groups):
      - levels: like compute_levels()'s result, EXCEPT that a node
        belonging to a fission group is given level len(loop_vars) --
        i.e. "only reachable at the very bottom of the main nest" is still
        literally true of where it's *referenced* from (see FissionGroup's
        docstring on why: even though a fission group's own auxiliary loop
        sits much shallower, at `depth`, a consumer outside the group can
        only safely load the scratch result back once every gap variable
        is genuinely back in scope in the main chain -- which, for the
        common two-dof-loop-plus-quadrature case this was written for, is
        only true at the innermost level anyway, since the gap variable is
        the innermost loop var itself). Every OTHER node's level is
        identical to compute_levels()'s -- this function only ever pulls
        nodes OUT of the "stuck at the deepest level" bucket compute_levels
        would otherwise put them in; it never changes the level of a node
        that already had a real prefix-covering depth.
      - fission_groups: one FissionGroup per distinct (depth, gap_vars)
        pairing that actually occurs, in the order their auxiliary loops
        should be emitted (by depth, and -- for groups sharing a depth --
        by which one is needed first in `nodes`' own dependency order, so a
        group whose nodes depend on another same-depth group's scratch
        output is never emitted before it).
    """
    loop_var_set = set(loop_vars)
    depends_on: dict[GraphNode, frozenset[str]] = {}
    graph = generate_graph(add_node.body)
    topo = list(graph.ordered_nodes(NodeOrder.leaves_first))
    for node in topo:
        deps = set(_direct_loop_deps(node, loop_var_set))
        for child in node.successors:
            deps |= depends_on[child]
        depends_on[node] = frozenset(deps)
    topo_index = {node: i for i, node in enumerate(topo)}

    # A node's "signature" is None if it's an ordinary prefix-hoistable node
    # (compute_levels' notion of level applies directly, no fission needed),
    # or (depth, gap_vars) if it needs fission.
    #
    # ArrayEntry is special-cased out of fission below. Fission (and the
    # alpha-equivalence scratch-sharing built on top of it, see
    # uflx_mlir.emit._alpha_signature) exists to avoid recomputing a
    # genuinely expensive shared subexpression (originally: geometry/
    # Jacobian terms, see this module's docstring). A bare ArrayEntry --
    # a table lookup with no further computation, e.g. an FE0 basis
    # value read -- has no such cost: reloading it at the point of use
    # costs exactly what reading a cached scratch value costs, so
    # fissioning one only adds a pointless copy-to-alloca loop. Verified
    # empirically (disassembly + timing A/B, correctness-checked via
    # np.testing.assert_allclose): excluding ArrayEntry from fission
    # dropped 250->222 instructions and gave a measured ~4% (1.0408x)
    # per-call speedup on a P3 stiffness kernel, with zero change to the
    # numeric result. It also matches FFCx's own generated C code, which
    # always re-reads its static table arrays directly at each point of
    # use rather than ever caching a plain table value.
    #
    # ArrayEntry nodes are pure leaves (no successors), so excluding them
    # here cannot change any other node's `deps`/gap computation -- only
    # this node's own signature/level assignment is affected, and forcing
    # levels[node] = len(loop_vars) below (innermost) mirrors what a
    # fission group would have assigned it anyway.
    signature: dict[GraphNode, tuple[int, tuple[str, ...]] | None] = {}
    levels: dict[GraphNode, int] = {}
    for node, deps in depends_on.items():
        depth = _prefix_depth(deps, loop_vars)
        gap_vars = tuple(v for v in loop_vars[depth:] if v in deps)
        if gap_vars and not isinstance(node, ArrayEntry):
            signature[node] = (depth, gap_vars)
        else:
            signature[node] = None
            levels[node] = len(loop_vars) if gap_vars else depth

    # Reverse of `successors` (a node's own operands) -- who *uses* each
    # node -- needed to tell whether a fission candidate's value ever
    # escapes its own group. A sentinel stands in for "used by
    # AddToLocalTensor itself" (add_node.body has no successors of its
    # own within this DAG), so a fissioned root is still correctly seen as
    # escaping its group.
    root_consumer = object()
    parents: dict[GraphNode, list[object]] = {node: [] for node in topo}
    for node in topo:
        for child in node.successors:
            parents[child].append(node)
    parents[add_node.body].append(root_consumer)

    groups: dict[tuple[int, tuple[str, ...]], list[GraphNode]] = {}
    for node in topo:
        sig = signature[node]
        if sig is not None:
            groups.setdefault(sig, []).append(node)

    fission_groups: list[FissionGroup] = []
    for (depth, gap_vars), nodes in groups.items():
        nodes_sorted = tuple(sorted(nodes, key=topo_index.__getitem__))
        scratch = tuple(
            node
            for node in nodes_sorted
            if any(
                parent is root_consumer
                or signature.get(cast(GraphNode, parent)) != (depth, gap_vars)
                for parent in parents[node]
            )
        )
        for node in nodes_sorted:
            levels[node] = len(loop_vars)
        fission_groups.append(FissionGroup(depth, gap_vars, nodes_sorted, scratch))

    fission_groups.sort(key=lambda g: (g.depth, min(topo_index[n] for n in g.nodes)))
    return levels, fission_groups


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
