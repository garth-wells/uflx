"""Generate MLIR for a UFLx form using the Python op-builder API.

This is deliberately NOT a from-scratch FEM-kernel generator: all of the
actual finite-element math (quadrature selection/tabulation, geometry
expansion, pull-back/push-forward, inner-product expansion) is done by
uflx_codegeneration's existing, generic, C-backend pipeline, reused
unchanged via `uflx_mlir.lowering.lower_form`. The only thing this module
adds is a replacement for the LAST step: instead of
uflx_codegeneration.c's monkey-patched `.generate_c()` string builder, it
walks the same fully-lowered graph and builds MLIR IR directly with
`mlir.ir.Operation.create`.

Two things keep the generated code from doing needless work:
  - memoization per node (a node shared by multiple parents is computed
    once).
  - explicit loop-level hoisting (see uflx_mlir.hoist): uflx_codegeneration
    nests the dof loops OUTSIDE the quadrature loop, which -- combined
    with naive emission that placed every node at the innermost loop
    unconditionally -- meant the WHOLE per-entry expression (including the
    cell's Jacobian/detJ/cofactor terms, invariant across every dof pair)
    was recomputed once per (dof, dof, quadrature-point) combination.
    generate_mlir_module() below instead (a) reorders the loop nest so the
    quadrature loop is outermost and (b) places each node of the
    AddToLocalTensor body at the shallowest loop depth its actual
    dependencies allow, computed directly from the graph rather than
    relying on the LLVM/MLIR optimizer to rediscover it after the fact --
    that was tried first and confirmed NOT to work (MLIR's generic
    loop-invariant-code-motion pass doesn't hoist memref.load-derived
    computation without stronger alias guarantees than are available
    here). See hoist.py's module docstring for the measured effect
    (~27x slower than FFCx per call before this, ~1.8x faster after, on a
    P3 tetrahedron stiffness kernel).

    Depth-based hoisting alone still has a blind spot: a node that depends
    on the quadrature point and one dof loop but not the OTHER dof loop
    sitting between them in the nesting (e.g. depends on {quadrature,
    trial-dof} but not test-dof, with loops nested quadrature -> test-dof
    -> trial-dof) can't be covered by any PREFIX of the loop nest, so
    hoist.compute_levels() alone would force it to the innermost level
    anyway -- costing O(ntest_dofs * ntrial_dofs * nquadrature) instead of
    the O(ntrial_dofs * nquadrature) its actual dependencies allow. This
    dominates at higher polynomial degree (larger ndofs), which is exactly
    what regressed a P6 tetrahedron stiffness kernel (84 dofs) to ~10%
    slower than FFCx despite the identical-shaped P3 kernel (20 dofs)
    being ~1.8x faster. hoist.compute_fission_plan() detects this "gap"
    pattern directly and reports it as a FissionGroup; _emit_fission_group
    below materializes each such group into its own small auxiliary loop
    -- covering just the skipped-over variable(s) -- storing results into
    a scratch memref and reading them back via a plain array lookup at
    the point they're actually needed, the same technique the
    hand-written kernels in mlir-kernels/python/generate_kernel.py already
    use (their numx_scratch/numy_scratch/numz_scratch buffers).

Deliberately built almost entirely on `Operation.create` plus the core
`mlir.ir` primitives (Context, Location, Module, InsertionPoint, Region,
Block, the Type/Attribute constructors) rather than the auto-generated
per-dialect Python wrapper classes (`arith.AddFOp`, `memref.LoadOp`,
`scf.ForOp`, `func.FuncOp`, ...). `Operation.create` only requires the
registered op name and its operand/attribute names, which don't drift
across LLVM releases the way the wrapper classes' auto-generated
constructor signatures occasionally do. scf.for's region/block and
func.func's entry block are therefore built by hand
(`Region.blocks.append(*arg_types)` then `InsertionPoint(that block)`),
which is more verbose than the `scf.ForOp`/`func.FuncOp` conveniences but
depends on nothing beyond the same core primitives.

Only covers what uflx_codegeneration's own pipeline currently covers: a
single cell (codim-0) integral, no coefficients/constants, a single
scalar element per function space (see `integrals_to_quadrature`'s own
asserts). Requires the `mlir` Python package (built from LLVM/MLIR source
with Python bindings enabled -- see this package's README.md); it is
deliberately not listed as a pip dependency since no pip wheel reliably
provides it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import basix
import numpy as np
from mlir.ir import (
    Context,
    DenseElementsAttr,
    F64Type,
    FlatSymbolRefAttr,
    FloatAttr,
    FunctionType,
    IndexType,
    InsertionPoint,
    IntegerAttr,
    Location,
    MemRefType,
    Module,
    Operation,
    RankedTensorType,
    StringAttr,
    TypeAttr,
    UnitAttr,
    Value,
)
from uflx.expressions import (
    Abs,
    Add,
    Div,
    Integer,
    Mult,
    Neg,
    RealScalar,
    Subtract,
)
from uflx.geometry import CoordinateDofComponent
from uflx.graphs import GraphNode
from uflx_codegeneration.nodes import AddToLocalTensor, ArrayEntry, FunctionCall, Loop, Variable
from uflx_codegeneration.quadrature import QuadratureLoop

from uflx_mlir.hoist import (
    FissionGroup,
    compute_fission_plan,
    reorder_quadrature_outermost,
    topo_order,
    walk_loop_chain,
)
from uflx_mlir.lowering import collect_int_constants, coordinate_shape, lower_form


@dataclass
class _OpCtx:
    """Mutable state threaded through op-builder emission."""

    a_shape: tuple[int, ...]
    coords_shape: tuple[int, int]
    table_shapes: dict[str, tuple[int, ...]]

    f64: Any
    index_t: Any

    index_const: dict[int, Value] = field(default_factory=dict)  # int -> index-typed Value
    global_val: dict[str, Value] = field(default_factory=dict)  # table name -> memref Value
    index_vars: dict[str, Value] = field(default_factory=dict)  # loop var name -> index Value
    zero_f64: Value | None = None

    # Fission-group scratch buffers (see hoist.compute_fission_plan). A node
    # is a key here only once its group's auxiliary loop has actually
    # stored its value -- see _emit_fission_group -- so that _emit_node
    # only takes the "load from scratch" path once that value genuinely
    # exists and is back in scope.
    scratch_group: dict[Any, FissionGroup] = field(default_factory=dict)
    scratch_buf: dict[Any, Value] = field(default_factory=dict)
    scratch_type: dict[Any, Any] = field(default_factory=dict)

    def resolve_index(self, idx: int | str) -> Value:
        """Resolve a Loop/ArrayEntry index entry to an `index`-typed SSA value.

        Args:
            idx: An int, a loop-variable name, or the literal string "0"
                uflx_codegeneration uses for a size-1 tensor axis.

        Returns:
            The corresponding `index`-typed Value.
        """
        if isinstance(idx, int):
            return self.index_const[idx]
        if isinstance(idx, str):
            if idx.lstrip("-").isdigit():
                return self.index_const[int(idx)]
            return self.index_vars[idx]
        raise NotImplementedError(f"Cannot resolve index of type {type(idx)}")


def _op1(name: str, result_type, operand) -> Value:
    return Operation.create(name, results=[result_type], operands=[operand]).results[0]


def _op2(name: str, result_type, lhs, rhs) -> Value:
    return Operation.create(name, results=[result_type], operands=[lhs, rhs]).results[0]


def _const_f64(ctx: _OpCtx, value: float) -> Value:
    return Operation.create(
        "arith.constant",
        results=[ctx.f64],
        attributes={"value": FloatAttr.get(ctx.f64, float(value))},
    ).results[0]


def _const_index(ctx: _OpCtx, value: int) -> Value:
    return Operation.create(
        "arith.constant",
        results=[ctx.index_t],
        attributes={"value": IntegerAttr.get(ctx.index_t, value)},
    ).results[0]


def _memref_load(memref_val: Value, indices: list[Value], elem_type) -> Value:
    return Operation.create(
        "memref.load", results=[elem_type], operands=[memref_val, *indices]
    ).results[0]


def _memref_store(value: Value, memref_val: Value, indices: list[Value]) -> None:
    Operation.create("memref.store", operands=[value, memref_val, *indices])


def _emit_node(node: GraphNode, cache: dict[Any, Value], ctx: _OpCtx) -> None:
    """Emit ops for one node of the AddToLocalTensor body's expression DAG.

    Every graph successor (child) of `node` already has an entry in
    `cache` -- guaranteed by driving this from hoist.topo_order()'s
    children-before-parents ordering, gated by hoist.compute_fission_plan()'s
    depth assignment (see generate_mlir_module below). Unlike a naive
    recursive walk, this never recurses into children -- they're looked
    up directly, since they're guaranteed already emitted. Ops are
    created at whatever the CURRENT InsertionPoint is, which the
    depth-driven driver in generate_mlir_module has already positioned at
    the right loop level.
    """
    if node in cache:
        return

    if node in ctx.scratch_group:
        # Computed ahead of time by a fission group's auxiliary loop (see
        # _emit_fission_group) -- read back here at its actual use site
        # with a fresh memref.load, not recomputed.
        group = ctx.scratch_group[node]
        idx = [ctx.resolve_index(v) for v in group.gap_vars]
        v = _memref_load(ctx.scratch_buf[node], idx, ctx.f64)
        cache[node] = v
        return

    if isinstance(node, (RealScalar, Integer)):
        v = _const_f64(ctx, node.value)
    elif isinstance(node, Neg):
        a = cache[node.argument]
        v = _op1("arith.negf", ctx.f64, a)
    elif isinstance(node, Abs):
        a = cache[node.argument]
        v = _op1("math.absf", ctx.f64, a)
    elif isinstance(node, Add):
        a = cache[node.first]
        b = cache[node.second]
        v = _op2("arith.addf", ctx.f64, a, b)
    elif isinstance(node, Subtract):
        a = cache[node.first]
        b = cache[node.second]
        v = _op2("arith.subf", ctx.f64, a, b)
    elif isinstance(node, Mult):
        a = cache[node.first]
        b = cache[node.second]
        v = _op2("arith.mulf", ctx.f64, a, b)
    elif isinstance(node, Div):
        a = cache[node.first]
        b = cache[node.second]
        v = _op2("arith.divf", ctx.f64, a, b)
    elif isinstance(node, CoordinateDofComponent):
        pt = ctx.resolve_index(node._point)
        comp = ctx.resolve_index(node._component)
        v = _memref_load(ctx.coords_val, [pt, comp], ctx.f64)  # type: ignore[attr-defined]
    elif isinstance(node, ArrayEntry):
        mem = ctx.global_val[node.array]
        idx = [ctx.resolve_index(i) for i in node.index]
        v = _memref_load(mem, idx, ctx.f64)
    elif isinstance(node, (Variable, FunctionCall)):
        raise NotImplementedError(
            f"{type(node).__name__} has no op-builder emitter -- this means the form "
            "needs a feature (coefficients/constants) this prototype doesn't cover yet."
        )
    else:
        raise NotImplementedError(f"No op-builder emitter for expression node type {type(node)}")

    cache[node] = v


def _chain_bounds(chain: list[tuple[GraphNode, str]], var: str) -> tuple[int, int]:
    """Return the (lo, hi) integer bounds of the loop bound to `var` in `chain`."""
    for loop_node, v in chain:
        if v != var:
            continue
        if isinstance(loop_node, QuadratureLoop):
            return 0, loop_node.rule.npoints
        if isinstance(loop_node, Loop):
            assert isinstance(loop_node.start, int) and isinstance(loop_node.end, int)
            assert loop_node.start == 0, (
                f"loop variable {var!r} starts at {loop_node.start} -- fission scratch buffers "
                "assume every loop starts at 0."
            )
            return loop_node.start, loop_node.end
        raise NotImplementedError(f"Unexpected loop node type {type(loop_node)} in chain")
    raise AssertionError(f"loop variable {var!r} not found in chain {chain!r}")


def _emit_fission_group(
    group: FissionGroup,
    chain: list[tuple[GraphNode, str]],
    cache: dict[Any, Value],
    ctx: _OpCtx,
) -> None:
    """Emit one fission group's auxiliary precompute loop (see hoist.FissionGroup).

    Builds a nested scf.for (one per group.gap_vars entry) at the CURRENT
    InsertionPoint, which _build_nest has already positioned at
    `group.depth` -- every loop_vars[:group.depth] is already open (so
    group.nodes' non-gap dependencies, already emitted at that shallower
    level, are safely in cache and in scope), and nothing group.nodes
    computes leaks out except through the scratch buffers.
    """
    shape = tuple(_chain_bounds(chain, v)[1] for v in group.gap_vars)
    mem_ty = MemRefType.get(list(shape), ctx.f64)
    for node in group.scratch:
        buf = Operation.create("memref.alloc", results=[mem_ty]).results[0]
        ctx.scratch_buf[node] = buf
        ctx.scratch_type[node] = mem_ty

    def build(remaining: tuple[str, ...]) -> None:
        if not remaining:
            for node in group.nodes:
                _emit_node(node, cache, ctx)
            for node in group.scratch:
                idx = [ctx.resolve_index(v) for v in group.gap_vars]
                _memref_store(cache[node], ctx.scratch_buf[node], idx)
            return
        var = remaining[0]
        lo, hi = _chain_bounds(chain, var)
        lb, ub, step = ctx.index_const[lo], ctx.index_const[hi], ctx.index_const[1]
        for_op = Operation.create("scf.for", operands=[lb, ub, step], regions=1)
        block = for_op.regions[0].blocks.append(ctx.index_t)
        with InsertionPoint(block):
            assert var not in ctx.index_vars, f"loop variable {var!r} shadows an outer one"
            ctx.index_vars[var] = block.arguments[0]
            build(remaining[1:])
            Operation.create("scf.yield")
        del ctx.index_vars[var]

    build(group.gap_vars)

    scratch_set = set(group.scratch)
    for node in group.nodes:
        # Every value just computed is scoped to the scf.for regions built
        # above and may not be referenced again by identity outside them.
        # A scratch node gets a fresh, properly scoped replacement the
        # next time it's needed (see _emit_node's scratch_group branch); a
        # pure local intermediate (not in group.scratch) is never
        # referenced again at all, since every one of its consumers is
        # inside this same group.
        del cache[node]
        if node in scratch_set:
            ctx.scratch_group[node] = group


def _build_nest(
    depth: int,
    chain: list[tuple[GraphNode, str]],
    levels: dict[GraphNode, int],
    topo: list[GraphNode],
    emitted: set[GraphNode],
    cache: dict[Any, Value],
    ctx: _OpCtx,
    add_node: AddToLocalTensor,
    groups_by_depth: dict[int, list[FissionGroup]],
) -> None:
    """Recursively build the (reordered, hoisted) loop nest around add_node's accumulation.

    At each depth, first emits (in dependency order) every not-yet-emitted
    body-DAG node whose level allows it here, then emits any fission group
    anchored at this depth (see hoist.FissionGroup / _emit_fission_group),
    then either recurses into the next scf.for (depth < len(chain)) or
    performs the final load/add/store into %A (depth == len(chain)).
    """
    for node in topo:
        if node in emitted:
            continue
        if levels[node] <= depth:
            _emit_node(node, cache, ctx)
            emitted.add(node)

    for group in groups_by_depth.get(depth, []):
        _emit_fission_group(group, chain, cache, ctx)
        scratch_set = set(group.scratch)
        for node in group.nodes:
            # Scratch nodes stay out of `emitted`: they're picked back up
            # (via _emit_node's scratch_group redirect) the next time the
            # topo scan above reaches them, at whatever depth they're
            # actually used. A pure local intermediate is done for good --
            # nothing outside this group ever references it -- so it must
            # be marked emitted, or a later depth's topo scan would try to
            # recompute it from since-deleted cache entries.
            if node not in scratch_set:
                emitted.add(node)

    if depth == len(chain):
        result = cache[add_node.body]
        idx = [ctx.resolve_index(i) for i in add_node.component]
        old = _memref_load(ctx.a_val, idx, ctx.f64)  # type: ignore[attr-defined]
        new = _op2("arith.addf", ctx.f64, old, result)
        _memref_store(new, ctx.a_val, idx)  # type: ignore[attr-defined]
        return

    loop_node, var = chain[depth]
    if isinstance(loop_node, QuadratureLoop):
        lo, hi = 0, loop_node.rule.npoints
    elif isinstance(loop_node, Loop):
        assert isinstance(loop_node.start, int) and isinstance(loop_node.end, int)
        lo, hi = loop_node.start, loop_node.end
    else:
        raise NotImplementedError(f"Unexpected loop node type {type(loop_node)} in chain")
    assert var not in ctx.index_vars, f"loop variable {var!r} shadows an outer one"

    lb = ctx.index_const[lo]
    ub = ctx.index_const[hi]
    step = ctx.index_const[1]
    for_op = Operation.create("scf.for", operands=[lb, ub, step], regions=1)
    block = for_op.regions[0].blocks.append(ctx.index_t)
    with InsertionPoint(block):
        ctx.index_vars[var] = block.arguments[0]
        _build_nest(depth + 1, chain, levels, topo, emitted, cache, ctx, add_node, groups_by_depth)
        Operation.create("scf.yield")
    del ctx.index_vars[var]


def _emit_zero_init(
    shape: tuple[int, ...], ctx: _OpCtx, indices: list[Value] | None = None
) -> None:
    indices = indices or []
    if len(indices) == len(shape):
        _memref_store(ctx.zero_f64, ctx.a_val, indices)  # type: ignore[attr-defined]
        return
    axis = len(indices)
    lb = ctx.index_const[0]
    ub = ctx.index_const[shape[axis]]
    step = ctx.index_const[1]
    for_op = Operation.create("scf.for", operands=[lb, ub, step], regions=1)
    block = for_op.regions[0].blocks.append(ctx.index_t)
    with InsertionPoint(block):
        _emit_zero_init(shape, ctx, [*indices, block.arguments[0]])
        Operation.create("scf.yield")


def generate_mlir_module(form, degree: int, kernel_name: str, cell: basix.CellType) -> Module:
    """Generate a self-contained MLIR module for a UFLx form's tabulate_tensor kernel.

    Built via the Python op-builder API rather than parsing MLIR text.

    Args:
        form: A UFLx form (as returned by e.g. `inner(grad(u), grad(v)) * dx`).
        degree: The polynomial degree used to size the quadrature rule.
        kernel_name: The symbol name to give the generated `func.func`.
        cell: The reference cell the form is integrated over.

    Returns:
        The built `mlir.ir.Module` (already verified with
        `module.operation.verify()`); its owning `Context` is kept alive
        via the Module's own reference to it.
    """
    ncoorddofs, tdim = coordinate_shape(form)
    tables, graph = lower_form(form, degree, cell)
    root = graph.root
    chain, add_node = walk_loop_chain(root)
    chain = reorder_quadrature_outermost(chain)
    loop_vars = [v for _, v in chain]
    a_shape = add_node.shape

    int_constants = collect_int_constants(root, a_shape)

    ctx_container = Context()
    with ctx_container, Location.unknown():
        module = Module.create()
        f64 = F64Type.get()
        index_t = IndexType.get()
        a_ty = MemRefType.get(list(a_shape), f64)
        coords_ty = MemRefType.get([ncoorddofs, tdim], f64)

        ctx = _OpCtx(
            a_shape=a_shape,
            coords_shape=(ncoorddofs, tdim),
            table_shapes={name: arr.shape for name, arr in tables.items()},
            f64=f64,
            index_t=index_t,
        )

        with InsertionPoint(module.body):
            table_types = {}
            for name in sorted(tables):
                arr = tables[name]
                ty = MemRefType.get(list(arr.shape), f64)
                table_types[name] = ty
                tensor_ty = RankedTensorType.get(list(arr.shape), f64)
                dense = DenseElementsAttr.get(
                    np.ascontiguousarray(arr, dtype=np.float64), type=tensor_ty
                )
                Operation.create(
                    "memref.global",
                    attributes={
                        "sym_name": StringAttr.get(name),
                        "sym_visibility": StringAttr.get("private"),
                        "type": TypeAttr.get(ty),
                        "initial_value": dense,
                        "constant": UnitAttr.get(),
                    },
                )

            func_ty = FunctionType.get([a_ty, coords_ty], [])
            func_op = Operation.create(
                "func.func",
                attributes={
                    "sym_name": StringAttr.get(kernel_name),
                    "function_type": TypeAttr.get(func_ty),
                    "llvm.emit_c_interface": UnitAttr.get(),
                },
                regions=1,
            )
            entry = func_op.regions[0].blocks.append(a_ty, coords_ty)
            with InsertionPoint(entry):
                ctx.a_val = entry.arguments[0]  # type: ignore[attr-defined]
                ctx.coords_val = entry.arguments[1]  # type: ignore[attr-defined]

                for i in sorted(int_constants):
                    ctx.index_const[i] = _const_index(ctx, i)
                ctx.zero_f64 = _const_f64(ctx, 0.0)

                for name in sorted(tables):
                    ctx.global_val[name] = Operation.create(
                        "memref.get_global",
                        results=[table_types[name]],
                        attributes={"name": FlatSymbolRefAttr.get(name)},
                    ).results[0]

                _emit_zero_init(a_shape, ctx)

                levels, fission_groups = compute_fission_plan(add_node, loop_vars)
                topo = topo_order(add_node)
                groups_by_depth: dict[int, list[FissionGroup]] = {}
                for group in fission_groups:
                    groups_by_depth.setdefault(group.depth, []).append(group)
                _build_nest(0, chain, levels, topo, set(), {}, ctx, add_node, groups_by_depth)

                Operation.create("func.return")

        module.operation.verify()

    return module
