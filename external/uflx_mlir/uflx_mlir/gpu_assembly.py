"""Experimental GPU-assembly entry-point kernel: one thread, one matrix entry.

Companion to emit.generate_mlir_module (which is left completely untouched
by this module -- see emit.py's thread_bindings/commit hooks, both
default-off). This targets the shape of
dolfinx-gpu-solvers/cpp/demo/mat_assembly/laplacian.h: no cells_per_block
batching (that is future work, deliberately deferred), one thread per
local (i, j) entry, geometry supplied precomputed rather than computed
from coordinates in-kernel (reusing uflx_mlir.geometry's existing affine
Poisson extraction -- this only ever succeeds for the same affine
tetrahedron Poisson/stiffness shape that already benefits from it on the
CPU path), and a CSR binary-search-and-atomic-add scatter matching
laplacian.h's own loop instruction for instruction (including its
Arowptr[row+1]-as-initial-maxidx convention) -- kept as plain, later-
refinable boilerplate rather than something smarter, by explicit request.

Deliberately NOT yet wrapped in an actual `gpu.func`/`gpu.launch`: the
(tx, ty, cell_id) triple this module's entry point takes are plain
`index`-typed function arguments, not `gpu.thread_id`/`gpu.block_id`
values. That makes the kernel body below callable and verifiable with an
ordinary CPU-only LLVM/MLIR build (no NVPTX/AMDGPU target needed) --
useful right now, since deriving (tx, ty, cell_id) from an actual
`gpu.launch` is a separate, smaller step that specifically needs the
GPU-target build to test against. Wrapping it is the natural next step
once that build is available.

Known, deliberate inefficiency: the fission mechanism that shares a
sub-expression across a whole DOF LOOP on the CPU path (see hoist.py)
still runs unchanged here, but its benefit assumes the axis it hoists
*past* is a real loop with many iterations to share across. Here both dof
axes are single thread-bound values, so a fission group's own auxiliary
loop degenerates to a trip-count-1 loop computing one scratch entry that
is immediately read back once -- correct, but pure overhead (an
oversized scratch buffer and a pointless store-then-load) compared to
just inlining the computation. Not fixed here: teaching
hoist.compute_fission_plan to recognize a thread-bound gap variable and
skip fissioning it entirely is a natural follow-up, kept out of this
change to keep the diff reviewable.

This module uses `scf.while`, `scf.if`, and `memref.atomic_rmw`
(kind `addf`) -- none of which anything else in this codebase has used
before (emit.py sticks to `scf.for`, `arith`, and
`memref.load/store/alloca/global`). An earlier draft used
`memref.generic_atomic_rmw`/`memref.atomic_yield` (a compare-and-swap
loop) instead, on the theory that it needed no risky "kind" attribute to
guess; running it against a real mlir build showed that reasoning was
backwards for a float accumulator: LLVM's `cmpxchg` (what
`generic_atomic_rmw` lowers to) only accepts integer or pointer operands,
so it failed to compile for an f64 memref ("operand #1 must be integer or
LLVM pointer type, but got 'f64'"). `memref.atomic_rmw`'s `addf` kind
instead lowers to a native `llvm.atomicrmw fadd`, which NVVM/PTX
implements as a real hardware atomic float add (`atom.add.f64`) -- both
correct and a better match for laplacian.h's own `atomicAdd` than a
hand-rolled CAS loop ever was. `scf.while`/`scf.if` region wiring and the
dynamic-memref-shape construction (`ShapedType.get_dynamic_size()`) were
the other parts flagged as least certain without a real build to test
against; both turned out fine once actually run -- see the module's own
test file for what could be checked without one (the binary-search
algorithm itself, independent of its MLIR encoding).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

import basix
import numpy as np
from mlir.dialects import gpu as gpu_d
from mlir.dialects import memref as memref_d
from mlir.ir import (
    Context,
    DenseElementsAttr,
    F64Type,
    FlatSymbolRefAttr,
    FunctionType,
    IndexType,
    InsertionPoint,
    IntegerAttr,
    IntegerType,
    Location,
    MemRefType,
    Module,
    Operation,
    RankedTensorType,
    ShapedType,
    StringAttr,
    SymbolRefAttr,
    TypeAttr,
    UnitAttr,
    UnrankedMemRefType,
    Value,
)
from uflx_codegeneration.nodes import Loop
from uflx_codegeneration.quadrature import QuadratureLoop

from uflx_mlir.emit import (
    _build_nest,
    _const_f64,
    _const_index,
    _emit_fission_scratch_allocas,
    _memref_load,
    _memref_store,
    _op1,
    _op2,
    _OpCtx,
)
from uflx_mlir.geometry import GeometryKernelSpec
from uflx_mlir.hoist import (
    FissionGroup,
    compute_fission_plan,
    reorder_quadrature_outermost,
    topo_order,
    walk_loop_chain,
)
from uflx_mlir.lowering import collect_int_constants, lower_form


@dataclass(frozen=True)
class CsrEntryLayout:
    """Shapes/strides the generated function's callers must respect.

    Attributes:
        ndofs: Local dofs per cell (both axes -- this first cut only
            covers equal test/trial spaces, matching laplacian.h).
        geometry_size: Length of the per-cell packed geometry array
            (`GeometryKernelSpec.output_size`).
        cell_dofs_stride: `cell_dofs` is a flat array, one cell's worth of
            `ndofs` entries at a time -- entry (cell, local_dof) lives at
            `cell * cell_dofs_stride + local_dof`. Equal to `ndofs`.
    """

    ndofs: int
    geometry_size: int
    cell_dofs_stride: int


def _cmpi(predicate: int, i1, lhs: Value, rhs: Value) -> Value:
    """Emit `arith.cmpi` with a raw predicate ordinal (generic-op form has no keyword).

    arith::CmpIPredicate ordinals (stable, long-standing): eq=0, ne=1,
    slt=2, sle=3, sgt=4, sge=5, ult=6, ule=7, ugt=8, uge=9.
    """
    return Operation.create(
        "arith.cmpi",
        results=[i1],
        operands=[lhs, rhs],
        attributes={"predicate": IntegerAttr.get(IntegerType.get_signless(64), predicate)},
    ).results[0]


def _binary_search_and_scatter(
    ctx: _OpCtx,
    value: Value,
    avals: Value,
    acols: Value,
    arowptr: Value,
    row_i32: Value,
    col_i32: Value,
    i32: Any,
) -> None:
    """Emit laplacian.h's binary-search-then-atomic-add CSR insertion, verbatim.

    Mirrors (instruction for instruction, including the
    Arowptr[row+1]-as-initial-maxidx convention):
        int minidx = Arowptr[row], maxidx = Arowptr[row + 1];
        while (minidx <= maxidx) {
            int idx = minidx + (maxidx - minidx) / 2;
            if (Acols[idx] < col) minidx = idx + 1;
            else if (Acols[idx] > col) maxidx = idx - 1;
            else { atomicAdd(&Avals[idx], value); break; }
        }
    The "break" is realised by forcing minidx=1, maxidx=0 on a match,
    which fails the loop's own `minidx <= maxidx` condition on the next
    check -- no separate "done" flag needed. Checked as plain Python
    (binary_search_reference.py) against this exact convention, including
    the transient one-past-row-end access maxidx starts at, before
    writing any of the below.
    """
    i1 = IntegerType.get_signless(1)
    one_i32 = Operation.create(
        "arith.constant", results=[i32], attributes={"value": IntegerAttr.get(i32, 1)}
    ).results[0]
    two_i32 = Operation.create(
        "arith.constant", results=[i32], attributes={"value": IntegerAttr.get(i32, 2)}
    ).results[0]
    zero_i32 = Operation.create(
        "arith.constant", results=[i32], attributes={"value": IntegerAttr.get(i32, 0)}
    ).results[0]

    row_idx = _op1("arith.index_cast", ctx.index_t, row_i32)
    row_idx_plus1 = _op2("arith.addi", ctx.index_t, row_idx, ctx.index_const[1])
    min0 = _memref_load(arowptr, [row_idx], i32)
    max0 = _memref_load(arowptr, [row_idx_plus1], i32)

    while_op = Operation.create("scf.while", results=[i32, i32], operands=[min0, max0], regions=2)

    before = while_op.regions[0].blocks.append(i32, i32)
    with InsertionPoint(before):
        cur_min, cur_max = before.arguments
        cond = _cmpi(3, i1, cur_min, cur_max)  # sle
        Operation.create("scf.condition", operands=[cond, cur_min, cur_max])

    after = while_op.regions[1].blocks.append(i32, i32)
    with InsertionPoint(after):
        cur_min, cur_max = after.arguments
        half = _op2("arith.divsi", i32, _op2("arith.subi", i32, cur_max, cur_min), two_i32)
        idx_i32 = _op2("arith.addi", i32, cur_min, half)
        idx = _op1("arith.index_cast", ctx.index_t, idx_i32)
        at_idx = _memref_load(acols, [idx], i32)
        lt = _cmpi(2, i1, at_idx, col_i32)  # slt
        gt = _cmpi(4, i1, at_idx, col_i32)  # sgt

        outer_if = Operation.create("scf.if", results=[i32, i32], operands=[lt], regions=2)
        then_blk = outer_if.regions[0].blocks.append()
        with InsertionPoint(then_blk):
            Operation.create(
                "scf.yield", operands=[_op2("arith.addi", i32, idx_i32, one_i32), cur_max]
            )
        else_blk = outer_if.regions[1].blocks.append()
        with InsertionPoint(else_blk):
            inner_if = Operation.create("scf.if", results=[i32, i32], operands=[gt], regions=2)
            inner_then = inner_if.regions[0].blocks.append()
            with InsertionPoint(inner_then):
                Operation.create(
                    "scf.yield", operands=[cur_min, _op2("arith.subi", i32, idx_i32, one_i32)]
                )
            inner_else = inner_if.regions[1].blocks.append()
            with InsertionPoint(inner_else):
                # AtomicRMWKind (arith::AtomicRMWKind, stable, long-standing
                # ordinals): addf=0, addi=1, assign=2, maximumf=3, maxs=4,
                # maxu=5, minimumf=6, mins=7, minu=8, mulf=9, muli=10, ori=11,
                # andi=12, maxnumf=13, minnumf=14. Built as a raw i64
                # IntegerAttr (confirmed against the real
                # AtomicRMWKindAttr registration, not guessed) rather than
                # importing the arith dialect's generated enum just for one
                # value.
                Operation.create(
                    "memref.atomic_rmw",
                    results=[ctx.f64],
                    operands=[value, avals, idx],
                    attributes={"kind": IntegerAttr.get(IntegerType.get_signless(64), 0)},
                )
                # Force minidx > maxidx so the while's own condition check
                # fails next time round -- this is the C loop's "break".
                Operation.create("scf.yield", operands=[one_i32, zero_i32])
            Operation.create("scf.yield", operands=list(inner_if.results))
        Operation.create("scf.yield", operands=list(outer_if.results))


def generate_csr_entry_module(
    form, degree: int, kernel_name: str, cell: basix.CellType
) -> tuple[Module, CsrEntryLayout]:
    """Generate a self-contained MLIR module for one (cell, i, j) CSR-entry kernel.

    See the module docstring: no batching yet, (cell_id, tx, ty) are plain
    `index` arguments (not derived from `gpu.*` ops), geometry is a
    caller-supplied packed array (not computed from coordinates), and the
    CSR insertion mirrors laplacian.h's binary search exactly.

    Only supports forms uflx_mlir.geometry.extract_affine_poisson_geometry
    can extract from (raises NotImplementedError otherwise -- this is a
    deliberate, documented scope limit, not a bug), and only equal
    test/trial local dof counts (also matching laplacian.h).

    Args:
        form: A UFLx form (as returned by e.g. `inner(grad(u), grad(v)) * dx`).
        degree: The polynomial degree used to size the quadrature rule.
        kernel_name: The symbol name to give the generated `func.func`.
        cell: The reference cell the form is integrated over.

    Returns:
        A tuple (module, layout): `module` is the built, verified
        `mlir.ir.Module`; `layout` records the shapes/strides callers must
        respect (see CsrEntryLayout).

    Raises:
        NotImplementedError: If the form's geometry can't be extracted, if
            the loop chain isn't the standard quadrature-plus-two-dof-axis
            shape, or if the two dof axes have unequal length.
    """
    tables, graph, geometry = lower_form(form, degree, cell)
    if geometry is None:
        raise NotImplementedError(
            "generate_csr_entry_module only supports forms whose geometry "
            "uflx_mlir.geometry.extract_affine_poisson_geometry can extract "
            "(currently: affine-tetrahedron Poisson/stiffness forms) -- see "
            "that module for why other forms are left alone."
        )
    assert isinstance(geometry, GeometryKernelSpec)

    root = graph.root
    chain, add_node = walk_loop_chain(root)
    chain = reorder_quadrature_outermost(chain)
    if len(chain) != 3 or not isinstance(chain[0][0], QuadratureLoop):
        raise NotImplementedError(
            "generate_csr_entry_module only supports the standard "
            "quadrature + two-dof-axis shape (a bilinear form with one "
            "test and one trial function loop)."
        )
    (_, quad_var), (dof_a_node, dof_a_var), (dof_b_node, dof_b_var) = chain
    assert isinstance(dof_a_node, Loop) and isinstance(dof_b_node, Loop)
    ndofs_a, ndofs_b = dof_a_node.end, dof_b_node.end
    # Loop.end is typed `int | str` (symbolic bounds are used elsewhere in
    # this codebase); narrow it the same way lowering.py's
    # collect_int_constants does before relying on it as a real int.
    if not isinstance(ndofs_a, int) or not isinstance(ndofs_b, int):
        raise NotImplementedError(
            "generate_csr_entry_module only supports forms with constant "
            f"(non-symbolic) local dof counts, got {ndofs_a!r} and {ndofs_b!r}."
        )
    if ndofs_a != ndofs_b:
        raise NotImplementedError(
            "generate_csr_entry_module only supports equal test/trial local "
            f"dof counts (got {ndofs_a} and {ndofs_b}), matching laplacian.h."
        )
    ndofs = ndofs_a
    loop_vars = [quad_var, dof_a_var, dof_b_var]

    int_constants = collect_int_constants(root, add_node.shape)
    int_constants.update(range(geometry.output_size))
    int_constants.add(2)

    ctx_container = Context()
    with ctx_container, Location.unknown():
        module = Module.create()
        f64 = F64Type.get()
        index_t = IndexType.get()
        i32 = IntegerType.get_signless(32)
        dyn = ShapedType.get_dynamic_size()

        avals_ty = MemRefType.get([dyn], f64)
        acols_ty = MemRefType.get([dyn], i32)
        arowptr_ty = MemRefType.get([dyn], i32)
        geometry_ty = MemRefType.get([geometry.output_size], f64)
        cell_dofs_ty = MemRefType.get([dyn], i32)

        ctx = _OpCtx(
            a_shape=add_node.shape,
            coords_shape=(0, 0),  # unused here: this kernel never touches coordinates
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

            func_ty = FunctionType.get(
                [
                    avals_ty,
                    acols_ty,
                    arowptr_ty,
                    geometry_ty,
                    cell_dofs_ty,
                    index_t,
                    index_t,
                    index_t,
                ],
                [],
            )
            func_op = Operation.create(
                "func.func",
                attributes={
                    "sym_name": StringAttr.get(kernel_name),
                    "function_type": TypeAttr.get(func_ty),
                    "llvm.emit_c_interface": UnitAttr.get(),
                },
                regions=1,
            )
            entry = func_op.regions[0].blocks.append(
                avals_ty,
                acols_ty,
                arowptr_ty,
                geometry_ty,
                cell_dofs_ty,
                index_t,
                index_t,
                index_t,
            )
            with InsertionPoint(entry):
                (
                    avals,
                    acols,
                    arowptr,
                    geometry_arg,
                    cell_dofs,
                    cell_id,
                    tx,
                    ty,
                ) = entry.arguments

                for i in sorted(int_constants):
                    ctx.index_const[i] = _const_index(ctx, i)
                ctx.zero_f64 = _const_f64(ctx, 0.0)

                ctx.geometry_val = geometry_arg
                ctx.thread_bindings = {dof_a_var: tx, dof_b_var: ty}

                for name in sorted(tables):
                    ctx.global_val[name] = Operation.create(
                        "memref.get_global",
                        results=[table_types[name]],
                        attributes={"name": FlatSymbolRefAttr.get(name)},
                    ).results[0]

                # Thread-private accumulator: sums this (cell, tx, ty)
                # entry's contribution across quadrature points via plain
                # (non-atomic -- nothing else touches this buffer)
                # load/add/store, the same idiom generate_mlir_module uses
                # for ctx.a_val, so only ONE atomic scatter happens total
                # rather than one per quadrature point.
                accum_ty = MemRefType.get([], f64)
                accum = Operation.create("memref.alloca", results=[accum_ty]).results[0]
                _memref_store(ctx.zero_f64, accum, [])

                def _accumulate(c: _OpCtx, result: Value) -> None:
                    old = _memref_load(accum, [], c.f64)
                    new = _op2("arith.addf", c.f64, old, result)
                    _memref_store(new, accum, [])

                ctx.commit = _accumulate

                levels, fission_groups = compute_fission_plan(add_node, loop_vars)
                _emit_fission_scratch_allocas(fission_groups, chain, ctx)

                topo = topo_order(add_node)
                groups_by_depth: dict[int, list[FissionGroup]] = {}
                for group in fission_groups:
                    groups_by_depth.setdefault(group.depth, []).append(group)
                _build_nest(0, chain, levels, topo, set(), {}, ctx, add_node, groups_by_depth)

                final = _memref_load(accum, [], f64)

                ndofs_const = ctx.index_const[ndofs]
                row_local = _op2(
                    "arith.addi", index_t, _op2("arith.muli", index_t, cell_id, ndofs_const), tx
                )
                col_local = _op2(
                    "arith.addi", index_t, _op2("arith.muli", index_t, cell_id, ndofs_const), ty
                )
                row_i32 = _memref_load(cell_dofs, [row_local], i32)
                col_i32 = _memref_load(cell_dofs, [col_local], i32)

                _binary_search_and_scatter(ctx, final, avals, acols, arowptr, row_i32, col_i32, i32)

                Operation.create("func.return")

        module.operation.verify()

    layout = CsrEntryLayout(ndofs=ndofs, geometry_size=geometry.output_size, cell_dofs_stride=ndofs)
    return module, layout


def gpu_module_name(kernel_name: str) -> str:
    """Return the gpu.module symbol a GPU-wrapped kernel is nested inside."""
    return f"{kernel_name}_gpu_module"


def gpu_launch_name(kernel_name: str) -> str:
    """Return the host-side gpu.launch_func wrapper's function symbol."""
    return f"{kernel_name}_launch"


def generate_csr_entry_gpu_module(
    form, degree: int, kernel_name: str, cell: basix.CellType
) -> tuple[Module, CsrEntryLayout]:
    """Generate generate_csr_entry_module's kernel wrapped in an actual
    gpu.func/gpu.module, with (tx, ty) derived from gpu.thread_id instead
    of taken as plain function arguments, plus a host-side
    gpu.launch_func wrapper -- the natural next step now that the
    NVPTX/AMDGPU-enabled LLVM build is available to test against.

    Scope, deliberately narrower than generate_csr_entry_module's own
    (already narrow) scope: one gpu.launch_func call handles exactly ONE
    cell (blocks in (1, 1, 1), threads in (ndofs, ndofs, 1)) -- cell_id is
    still a plain kernel argument, not derived from gpu.block_id. Batching
    many cells into a single launch (gridDim.x = ncells, cell_id =
    gpu.block_id x) is future work: it would also require geometry to
    become a full per-cell-indexed array rather than the single
    6-element slice used here (and here reused completely unchanged from
    generate_csr_entry_module), which is out of scope for this change --
    see this module's other docstrings for the same "hold off on
    batching for now" decision applied to the non-GPU-wrapped kernel.

    This is new, untested against a real mlir build with GPU-dialect
    support: gpu.func/gpu.module/gpu.launch_func construction (including
    the sym_name/gpu.kernel attributes, set directly via
    Operation.attributes[...] rather than through generated builder
    keyword arguments, since GPUFuncOp/GPUModuleOp's Python builders
    don't expose sym_name as a constructor parameter -- confirmed against
    this build's actual generated _gpu_ops_gen.py, not guessed) has not
    been exercised before in this codebase. module.operation.verify() is
    called before returning, same as generate_csr_entry_module, so a
    construction mistake should surface immediately rather than silently
    misbuilding. Lowering this further to real NVVM/ROCDL
    (convert-gpu-to-nvvm / convert-gpu-to-rocdl, then
    gpu-module-to-binary) is a separate step not attempted here -- this
    only builds and verifies the gpu-dialect IR itself.

    The host launch wrapper also brackets its gpu.launch_func call with
    gpu.host_register/gpu.host_unregister on every memref argument.
    Without this, real execution (as opposed to just compiling to PTX/
    SASS -- see lower_module_to_nvvm) would be broken: gpu.launch_func's
    own lowering hands kernelOperands' raw pointers straight to
    mgpuLaunchKernel with no device allocation or host->device copy of
    its own, so an ordinary host-heap buffer would be an invalid pointer
    once dereferenced by device code. gpu.host_register maps existing
    host memory into the device address space in place instead (see its
    op doc in GPUOps.td) -- confirmed against GPUToLLVMConversion.cpp's
    mgpuMemHostRegisterMemRef call, not guessed. It only accepts
    unranked memrefs, hence the memref.cast before each call.

    Args/Returns/Raises: see generate_csr_entry_module -- identical, plus
    the returned module additionally contains the gpu.module (named
    gpu_module_name(kernel_name)) and the host launch wrapper (named
    gpu_launch_name(kernel_name)), both taking the same
    (Avals, Acols, Arowptr, geometry, cell_dofs, cell_id) argument list.
    """
    tables, graph, geometry = lower_form(form, degree, cell)
    if geometry is None:
        raise NotImplementedError(
            "generate_csr_entry_gpu_module only supports forms whose geometry "
            "uflx_mlir.geometry.extract_affine_poisson_geometry can extract "
            "(currently: affine-tetrahedron Poisson/stiffness forms) -- see "
            "that module for why other forms are left alone."
        )
    assert isinstance(geometry, GeometryKernelSpec)

    root = graph.root
    chain, add_node = walk_loop_chain(root)
    chain = reorder_quadrature_outermost(chain)
    if len(chain) != 3 or not isinstance(chain[0][0], QuadratureLoop):
        raise NotImplementedError(
            "generate_csr_entry_gpu_module only supports the standard "
            "quadrature + two-dof-axis shape (a bilinear form with one "
            "test and one trial function loop)."
        )
    (_, quad_var), (dof_a_node, dof_a_var), (dof_b_node, dof_b_var) = chain
    assert isinstance(dof_a_node, Loop) and isinstance(dof_b_node, Loop)
    ndofs_a, ndofs_b = dof_a_node.end, dof_b_node.end
    # Loop.end is typed `int | str`; narrow it before relying on it as a
    # real int (see the matching check in generate_csr_entry_module).
    if not isinstance(ndofs_a, int) or not isinstance(ndofs_b, int):
        raise NotImplementedError(
            "generate_csr_entry_gpu_module only supports forms with constant "
            f"(non-symbolic) local dof counts, got {ndofs_a!r} and {ndofs_b!r}."
        )
    if ndofs_a != ndofs_b:
        raise NotImplementedError(
            "generate_csr_entry_gpu_module only supports equal test/trial "
            f"local dof counts (got {ndofs_a} and {ndofs_b}), matching laplacian.h."
        )
    ndofs = ndofs_a
    loop_vars = [quad_var, dof_a_var, dof_b_var]

    int_constants = collect_int_constants(root, add_node.shape)
    int_constants.update(range(geometry.output_size))
    int_constants.add(2)

    ctx_container = Context()
    with ctx_container, Location.unknown():
        module = Module.create()
        # gpu.launch_func requires its closest surrounding module to carry
        # this attribute (verifier: "expected the closest surrounding
        # module to have the 'gpu.container_module' attribute") -- caught
        # by running this against a real mlir build.
        module.operation.attributes["gpu.container_module"] = UnitAttr.get()
        f64 = F64Type.get()
        index_t = IndexType.get()
        i32 = IntegerType.get_signless(32)
        dyn = ShapedType.get_dynamic_size()

        avals_ty = MemRefType.get([dyn], f64)
        acols_ty = MemRefType.get([dyn], i32)
        arowptr_ty = MemRefType.get([dyn], i32)
        geometry_ty = MemRefType.get([geometry.output_size], f64)
        cell_dofs_ty = MemRefType.get([dyn], i32)
        kernel_arg_types = [avals_ty, acols_ty, arowptr_ty, geometry_ty, cell_dofs_ty, index_t]

        ctx = _OpCtx(
            a_shape=add_node.shape,
            coords_shape=(0, 0),  # unused here: this kernel never touches coordinates
            table_shapes={name: arr.shape for name, arr in tables.items()},
            f64=f64,
            index_t=index_t,
        )

        with InsertionPoint(module.body):
            gmod_name = gpu_module_name(kernel_name)
            # GPUModuleOp's generated Python constructor has changed
            # across LLVM checkouts of the same release/18.x branch: an
            # older one (what this was originally written/verified
            # against) took no sym_name parameter at all -- set via
            # .attributes[...] afterward, same as GPUFuncOp below -- a
            # newer one requires it as a constructor argument outright
            # ("missing 1 required positional argument: 'sym_name'",
            # caught by running this against eng-nvidia's rebuilt LLVM).
            # Handle both rather than pin to whichever this happens to
            # be built against.
            try:
                gpu_module_op = gpu_d.GPUModuleOp(sym_name=StringAttr.get(gmod_name))
            except TypeError:
                gpu_module_op = gpu_d.GPUModuleOp()
                gpu_module_op.attributes["sym_name"] = StringAttr.get(gmod_name)
            gpu_body = gpu_module_op.regions[0].blocks.append()
            with InsertionPoint(gpu_body):
                # gpu.module is IsolatedFromAbove and its own SymbolTable
                # (GPUOps.td's GPU_GPUModuleOp traits): a kernel compiled
                # from it must be self-contained, so the FE/quadrature
                # table globals a memref.get_global inside the kernel body
                # resolves have to live IN HERE too, as siblings of
                # gpu.func -- unlike generate_csr_entry_module's plain
                # func.func, which can declare them at the outer module.
                # (Caught by running this against a real mlir build:
                # 'memref.get_global' op 'FE2' does not reference a valid
                # global memref.)
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

                kernel_func_ty = FunctionType.get(kernel_arg_types, [])
                gpu_func_op = gpu_d.GPUFuncOp(TypeAttr.get(kernel_func_ty))
                gpu_func_op.attributes["sym_name"] = StringAttr.get(kernel_name)
                # See GPUBase.td's getKernelFuncAttrName(): a gpu.func's
                # "kernel" keyword in textual IR is exactly this unit
                # attribute, not a constructor parameter.
                gpu_func_op.attributes["gpu.kernel"] = UnitAttr.get()
                entry = gpu_func_op.regions[0].blocks.append(*kernel_arg_types)
                with InsertionPoint(entry):
                    avals, acols, arowptr, geometry_arg, cell_dofs, cell_id = entry.arguments

                    for i in sorted(int_constants):
                        ctx.index_const[i] = _const_index(ctx, i)
                    ctx.zero_f64 = _const_f64(ctx, 0.0)

                    ctx.geometry_val = geometry_arg
                    tx = gpu_d.thread_id(gpu_d.Dimension.x)
                    ty = gpu_d.thread_id(gpu_d.Dimension.y)
                    ctx.thread_bindings = {dof_a_var: tx, dof_b_var: ty}

                    for name in sorted(tables):
                        ctx.global_val[name] = Operation.create(
                            "memref.get_global",
                            results=[table_types[name]],
                            attributes={"name": FlatSymbolRefAttr.get(name)},
                        ).results[0]

                    # Same thread-private accumulator idiom as
                    # generate_csr_entry_module -- see there for why.
                    accum_ty = MemRefType.get([], f64)
                    accum = Operation.create("memref.alloca", results=[accum_ty]).results[0]
                    _memref_store(ctx.zero_f64, accum, [])

                    def _accumulate(c: _OpCtx, result: Value) -> None:
                        old = _memref_load(accum, [], c.f64)
                        new = _op2("arith.addf", c.f64, old, result)
                        _memref_store(new, accum, [])

                    ctx.commit = _accumulate

                    levels, fission_groups = compute_fission_plan(add_node, loop_vars)
                    _emit_fission_scratch_allocas(fission_groups, chain, ctx)

                    topo = topo_order(add_node)
                    groups_by_depth: dict[int, list[FissionGroup]] = {}
                    for group in fission_groups:
                        groups_by_depth.setdefault(group.depth, []).append(group)
                    _build_nest(0, chain, levels, topo, set(), {}, ctx, add_node, groups_by_depth)

                    final = _memref_load(accum, [], f64)

                    ndofs_const = ctx.index_const[ndofs]
                    row_local = _op2(
                        "arith.addi", index_t, _op2("arith.muli", index_t, cell_id, ndofs_const), tx
                    )
                    col_local = _op2(
                        "arith.addi", index_t, _op2("arith.muli", index_t, cell_id, ndofs_const), ty
                    )
                    row_i32 = _memref_load(cell_dofs, [row_local], i32)
                    col_i32 = _memref_load(cell_dofs, [col_local], i32)

                    _binary_search_and_scatter(
                        ctx, final, avals, acols, arowptr, row_i32, col_i32, i32
                    )

                    gpu_d.ReturnOp(operands_=[])
                gpu_d.ModuleEndOp()

            # Host-side launch wrapper: one gpu.launch_func call, one cell
            # per launch (see this function's docstring for why).
            launch_name = gpu_launch_name(kernel_name)
            launch_func_ty = FunctionType.get(kernel_arg_types, [])
            launch_op = Operation.create(
                "func.func",
                attributes={
                    "sym_name": StringAttr.get(launch_name),
                    "function_type": TypeAttr.get(launch_func_ty),
                    "llvm.emit_c_interface": UnitAttr.get(),
                },
                regions=1,
            )
            launch_entry = launch_op.regions[0].blocks.append(*kernel_arg_types)
            with InsertionPoint(launch_entry):
                (
                    h_avals,
                    h_acols,
                    h_arowptr,
                    h_geometry,
                    h_cell_dofs,
                    h_cell_id,
                ) = launch_entry.arguments

                # gpu.launch_func's own lowering (GPUToLLVMConversion.cpp's
                # LaunchFuncOpLowering) does nothing more than hand the
                # kernelOperands' raw pointers straight to mgpuLaunchKernel
                # -- it does NOT allocate device memory or copy host data
                # there first. Without registering these host buffers, the
                # kernel would dereference ordinary host malloc'd pointers
                # from device code: undefined behaviour (illegal memory
                # access) on a discrete GPU. gpu.host_register (lowering to
                # mgpuMemHostRegisterMemRef, confirmed against
                # GPUToLLVMConversion.cpp) maps existing host memory into
                # the device address space in place instead -- no separate
                # device allocation or explicit copy needed, and per its
                # own op doc host writes made before this point (all of
                # ours are, since the caller fills these buffers before
                # invoking this function) are guaranteed visible to a
                # kernel launched afterwards, with device writes visible
                # back on the host once the (synchronous, no async token)
                # gpu.launch_func below returns. It only accepts unranked
                # memrefs (GPUOps.td's GPU_HostRegisterOp: AnyUnrankedMemRef)
                # so each ranked argument is memref.cast to that first.
                unranked_f64 = UnrankedMemRefType.get(f64, None)
                unranked_i32 = UnrankedMemRefType.get(i32, None)
                for value, unranked_ty in (
                    (h_avals, unranked_f64),
                    (h_acols, unranked_i32),
                    (h_arowptr, unranked_i32),
                    (h_geometry, unranked_f64),
                    (h_cell_dofs, unranked_i32),
                ):
                    gpu_d.HostRegisterOp(memref_d.CastOp(unranked_ty, value).result)

                # Fresh constants scoped to this block -- ctx.index_const's
                # entries were built inside the gpu.func above and are not
                # valid SSA values here.
                one = Operation.create(
                    "arith.constant",
                    results=[index_t],
                    attributes={"value": IntegerAttr.get(index_t, 1)},
                ).results[0]
                ndofs_const_host = Operation.create(
                    "arith.constant",
                    results=[index_t],
                    attributes={"value": IntegerAttr.get(index_t, ndofs)},
                ).results[0]

                gpu_d.LaunchFuncOp(
                    asyncToken=None,
                    asyncDependencies=[],
                    kernel=SymbolRefAttr.get([gmod_name, kernel_name]),
                    gridSizeX=one,
                    gridSizeY=one,
                    gridSizeZ=one,
                    blockSizeX=ndofs_const_host,
                    blockSizeY=ndofs_const_host,
                    blockSizeZ=one,
                    kernelOperands=[h_avals, h_acols, h_arowptr, h_geometry, h_cell_dofs, h_cell_id],
                )

                for value, unranked_ty in (
                    (h_avals, unranked_f64),
                    (h_acols, unranked_i32),
                    (h_arowptr, unranked_i32),
                    (h_geometry, unranked_f64),
                    (h_cell_dofs, unranked_i32),
                ):
                    gpu_d.HostUnregisterOp(memref_d.CastOp(unranked_ty, value).result)

                Operation.create("func.return")

        module.operation.verify()

    layout = CsrEntryLayout(ndofs=ndofs, geometry_size=geometry.output_size, cell_dofs_stride=ndofs)
    return module, layout


def lower_module_to_nvvm(
    module: Module,
    *,
    cubin_chip: str = "sm_80",
    cubin_format: str = "isa",
) -> None:
    """Lower a gpu.module-containing module (see generate_csr_entry_gpu_module)
    all the way to compiled NVVM/PTX, in place, via MLIR's own bundled
    `gpu-lower-to-nvvm-pipeline` (mlir/lib/Dialect/GPU/Pipelines/
    GPUToNVVMPipeline.cpp) rather than hand-assembling the individual
    conversion passes.

    That pipeline is only compiled in when the NVPTX backend is part of
    the build (MLIR_ENABLE_CUDA_CONVERSIONS, auto-derived in mlir's own
    CMakeLists.txt from "NVPTX" being in LLVM_TARGETS_TO_BUILD --
    confirmed true for the mlir-kernels build this was written against),
    and is registered automatically the moment `mlir.ir` is imported
    (mlirRegisterAllPasses() runs on module load -- see
    lib/Bindings/Python/RegisterEverything.cpp), so no separate
    registration call is needed here.

    cubin_format defaults to "isa" (plain PTX assembly text) rather than
    the pipeline's own default "fatbin": producing an actual cubin/fatbin
    requires shelling out to NVIDIA's `ptxas` (part of the CUDA Toolkit,
    found via $PATH or NVVM::getCUDAToolkitPath()'s env vars), which a
    plain development machine -- e.g. the Mac this was developed on --
    won't have installed; PTX text only needs LLVM's own compiled-in
    NVPTX backend, already confirmed present. Switch to "fatbin"/"bin" on
    a machine with the CUDA Toolkit installed to get an actually loadable
    binary rather than assembly text (see ModuleToBinary.cpp's
    CompilationTarget::{Assembly,Binary,Fatbin,Offload} for the accepted
    values -- both "isa"/"assembly" and "bin"/"binary"/"fatbin"/"fatbinary"
    spellings work).

    This exercises far more of the toolchain than anything else in this
    module: real NVPTX backend codegen, not just IR construction and
    module.operation.verify(). Running it against a real build (with
    cubin_format="fatbin", the original default) got all the way through
    every conversion pass -- including the arith/memref-to-llvm lowering
    of this module's memref.atomic_rmw<addf> into a native
    `llvm.atomicrmw fadd` -- and only failed at the very last step
    because `ptxas` wasn't on that machine; switching to "isa" avoids
    that step entirely. cubin_chip defaults to "sm_80" (Ampere) rather
    than the pipeline's own stale "sm_50" default; override it to match
    whatever GPU the compiled output is actually meant to target.

    No ROCDL/AMDGPU equivalent exists yet: unlike NVVM, this LLVM
    checkout has no bundled "gpu-lower-to-rocdl-pipeline" convenience
    pass (mlir/lib/Dialect/GPU/Pipelines/ contains only
    GPUToNVVMPipeline.cpp) -- an AMDGPU path would need the individual
    convert-gpu-to-rocdl / gpu-rocdl-attach-target / gpu-module-to-binary
    passes assembled by hand, not attempted here.

    Args:
        module: A module built by generate_csr_entry_gpu_module (or
            anything else containing a gpu.module targeting NVVM),
            modified in place.
        cubin_chip: The target NVPTX chip, e.g. "sm_70", "sm_80", "sm_90".
        cubin_format: The pipeline's `cubin-format` option -- "isa"
            (default here; human-readable PTX assembly, no `ptxas`
            needed), "fatbin" (the pipeline's own default) or "bin" (a
            single cubin, no fatbin wrapper) -- the latter two need
            `ptxas` on $PATH or findable via NVVM::getCUDAToolkitPath().
    """
    from mlir.passmanager import PassManager

    pipeline = (
        "builtin.module(gpu-lower-to-nvvm-pipeline{"
        f"cubin-chip={cubin_chip} cubin-format={cubin_format}"
        "})"
    )
    # generate_csr_entry_gpu_module's own `with ctx_container, ...:` block
    # (which made ctx_container the default/current context) has already
    # exited by the time a caller gets `module` back, so PassManager needs
    # the context passed explicitly here rather than relying on a
    # surrounding `with Context():` -- caught by running this against a
    # real mlir build ("An MLIR function requires a Context but none was
    # provided...").
    pm = PassManager.parse(pipeline, context=module.context)
    pm.run(module.operation)


def extract_ptx_text(module: Module) -> str:
    """Pull the raw PTX assembly text back out of a module already
    compiled by lower_module_to_nvvm(..., cubin_format="isa") (the
    default), so it can be handed to a real `ptxas` on a machine that
    actually has an NVIDIA GPU and the CUDA Toolkit -- this dev machine's
    build only got as far as PTX text (see lower_module_to_nvvm's own
    docstring for why).

    gpu-module-to-binary attaches its compiled output as a `#gpu.object<>`
    attribute inside a `gpu.binary` op it creates in place of the
    original `gpu.module`. There's no dedicated Python wrapper class for
    that attribute in this build (unlike most other GPU dialect
    attributes), so this recovers the string by hand instead: MLIR's
    string-attribute printer (llvm::printEscapedString, in
    llvm/lib/Support/StringExtras.cpp -- confirmed against that source,
    not guessed) prints every character that isn't printable-and-not-a-
    quote as a `\\XX` two-hex-digit escape -- including a literal `"`
    itself, which is why a plain `"([^"]*)"` regex is safe here: the
    payload can never contain an unescaped quote to confuse the string's
    boundary. Picks the LONGEST quoted string in the gpu.binary op's own
    printed text (rather than assuming a fixed position among the
    target/format/object fields) since the PTX payload is always far
    longer than the other string attributes next to it (e.g. the "isa"
    format tag or the "sm_80" chip name).

    Args:
        module: A module already processed by lower_module_to_nvvm.

    Returns:
        The decoded PTX assembly text.

    Raises:
        ValueError: If no gpu.binary op is found (module wasn't compiled
            yet, or compiled with a target other than NVVM) or it
            contains no string attribute to recover.
    """
    binary_op = None
    for op in module.body:
        if op.operation.name == "gpu.binary":
            binary_op = op
            break
    if binary_op is None:
        raise ValueError(
            "no gpu.binary op found in this module -- call "
            "lower_module_to_nvvm(module, ...) first"
        )

    text = str(binary_op.operation)
    candidates = re.findall(r'"([^"]*)"', text)
    if not candidates:
        raise ValueError("gpu.binary op has no quoted string attribute to recover")
    escaped = max(candidates, key=len)

    chars = []
    i = 0
    while i < len(escaped):
        c = escaped[i]
        if c == "\\":
            if escaped[i + 1] == "\\":
                chars.append("\\")
                i += 2
            else:
                chars.append(chr(int(escaped[i + 1 : i + 3], 16)))
                i += 3
        else:
            chars.append(c)
            i += 1
    # The serialized object string is null-terminated internally (an
    # artifact of how gpu-module-to-binary stores it, not part of the
    # actual PTX) -- confirmed by running this against a real build,
    # which decoded a trailing "\x00" every time. Strip it so callers
    # get exactly the PTX text, nothing else.
    return "".join(chars).rstrip("\x00")
