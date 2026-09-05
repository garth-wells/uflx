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

This module uses `scf.while`, `scf.if`, and
`memref.generic_atomic_rmw`/`memref.atomic_yield` -- none of which
anything else in this codebase has used before (emit.py sticks to
`scf.for`, `arith`, and `memref.load/store/alloca/global`). It was
written without access to a real `mlir.ir` build to verify against; the
`scf.while`/`scf.if` region wiring and the dynamic-memref-shape
construction (`ShapedType.get_dynamic_size()`) are the parts most likely
to need a fix once actually run -- see the module's own test file for
what could be checked without one (the binary-search algorithm itself,
independent of its MLIR encoding).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import basix
import numpy as np
from mlir.dialects import gpu as gpu_d
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
                rmw = Operation.create(
                    "memref.generic_atomic_rmw",
                    results=[ctx.f64],
                    operands=[avals, idx],
                    regions=1,
                )
                rmw_block = rmw.regions[0].blocks.append(ctx.f64)
                with InsertionPoint(rmw_block):
                    (old,) = rmw_block.arguments
                    new = _op2("arith.addf", ctx.f64, old, value)
                    Operation.create("memref.atomic_yield", operands=[new])
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
                Operation.create("func.return")

        module.operation.verify()

    layout = CsrEntryLayout(ndofs=ndofs, geometry_size=geometry.output_size, cell_dofs_stride=ndofs)
    return module, layout
