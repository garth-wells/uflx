"""Cooperative quadrature evaluation followed by test-DOF contraction on the GPU."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import cast

from mlir.dialects import gpu as gpu_d
from mlir.ir import (
    ArrayAttr,
    Attribute,
    DenseI32ArrayAttr,
    DenseI64ArrayAttr,
    DictAttr,
    FlatSymbolRefAttr,
    InsertionPoint,
    IntegerAttr,
    IntegerType,
    MemRefType,
    Operation,
    ShapedType,
    Value,
)
from uflx.graphs import GraphNode
from uflx_codegeneration.nodes import FunctionCall

from uflx_mlir.emit import (
    _const_f64,
    _const_index,
    _emit_coefficient_expr,
    _emit_node,
    _memref_load,
    _memref_store,
    _op1,
    _op2,
)
from uflx_mlir.gpu_assembly import _cmpi
from uflx_mlir.hoist import compute_levels, topo_order


def cooperative_memory_bytes(add_node, chain, tables, coefficient_size, cells_per_block):
    """Estimate static shared storage, including per-allocation alignment padding."""
    levels = compute_levels(add_node, [var for _, var in chain])
    frontier = {
        child
        for node in topo_order(add_node)
        if levels[node] == 2
        for child in node.successors
        if levels[child] < 2
    }
    unique = {(arr.shape, arr.dtype.str, arr.tobytes()): arr.nbytes for arr in tables.values()}
    nq = chain[0][0].rule.npoints
    return (
        sum(unique.values())
        + 15 * (len(unique) + 2)
        + 8 * cells_per_block * (coefficient_size + len(frontier) * nq)
    )


def emit_cooperative_action(
    kernel, arg_types, ctx, layout, chain, add_node, table_types, local_coeff, tables, functions
):
    """Share quadrature-only graph values across test DOFs, preserving the form."""
    loop_vars = [var for _, var in chain]
    levels = compute_levels(add_node, loop_vars)
    topo = topo_order(add_node)
    # These edges cross from quadrature-only evaluation to test-dependent work.
    frontier_set = {
        child
        for node in topo
        if levels[node] == 2
        for child in node.successors
        if levels[child] < 2
    }
    frontier = [node for node in topo if node in frontier_set]
    if not frontier:
        raise NotImplementedError("No quadrature-to-test frontier in cooperative action")
    space = Attribute.parse("#gpu.address_space<workgroup>")
    coeff_ty = MemRefType.get(
        [layout.coefficient_size, layout.cells_per_block], ctx.f64, memory_space=space
    )
    scratch_ty = MemRefType.get(
        [len(frontier), layout.quadrature_points, layout.cells_per_block],
        ctx.f64,
        memory_space=space,
    )
    unique = {}
    aliases = {}
    keys = {}
    for name, array in tables.items():
        key = (array.shape, array.dtype.str, array.tobytes())
        representative = keys.setdefault(key, name)
        aliases[name] = representative
        if representative == name:
            unique[name] = array
    shared_types = {
        name: MemRefType.get(list(arr.shape), ctx.f64, memory_space=space)
        for name, arr in unique.items()
    }
    shared_bytes = 8 * (
        layout.coefficient_size * layout.cells_per_block
        + len(frontier) * layout.quadrature_points * layout.cells_per_block
    )
    shared_bytes += sum(arr.nbytes for arr in unique.values())
    if shared_bytes > 48 * 1024:
        raise NotImplementedError("Cooperative action requires more than 48 KiB of shared memory")
    kernel.attributes["workgroup_attributions"] = IntegerAttr.get(
        IntegerType.get_signless(64), 2 + len(shared_types)
    )
    # Explicit alignment is required for vectorized shared double loads.
    kernel.attributes["workgroup_attrib_attrs"] = ArrayAttr.get(
        [
            DictAttr.get({"llvm.align": IntegerAttr.get(IntegerType.get_signless(64), 16)})
            for _ in range(2 + len(shared_types))
        ]
    )
    entry = kernel.regions[0].blocks.append(
        *arg_types, coeff_ty, scratch_ty, *shared_types.values()
    )
    with InsertionPoint(entry):
        args = cast(Sequence[Value], entry.arguments)
        output, metric, coefficients, dofmap, shared_coeffs, scratch = args[:6]
        shared_tables = dict(zip(shared_types, args[6:]))
        for value in sorted(
            {
                0,
                1,
                6,
                layout.ndofs,
                layout.quadrature_points,
                layout.coefficient_size,
                layout.cells_per_block,
                *range(len(frontier)),
                *range(6),
                *range(
                    0, layout.coefficient_size + layout.threads_per_cell, layout.threads_per_cell
                ),
            }
        ):
            ctx.index_const[value] = _const_index(ctx, value)
        ctx.zero_f64 = _const_f64(ctx, 0.0)
        c = ctx.index_const
        lane = cast(Value, gpu_d.thread_id(gpu_d.Dimension.y))
        cell_local = cast(Value, gpu_d.thread_id(gpu_d.Dimension.x))
        cell = _op2(
            "arith.addi",
            ctx.index_t,
            _op2(
                "arith.muli",
                ctx.index_t,
                gpu_d.block_id(gpu_d.Dimension.x),
                c[layout.cells_per_block],
            ),
            cell_local,
        )
        length = Operation.create(
            "memref.dim", results=[ctx.index_t], operands=[metric, c[0]]
        ).results[0]
        ncells = _op2("arith.divui", ctx.index_t, length, c[6])
        i1 = IntegerType.get_signless(1)
        valid_cell = _cmpi(6, i1, cell, ncells)
        for name, ty in table_types.items():
            ctx.global_val[name] = Operation.create(
                "memref.get_global", results=[ty], attributes={"name": FlatSymbolRefAttr.get(name)}
            ).results[0]
        block_threads = layout.threads_per_cell * layout.cells_per_block
        linear_thread = _op2(
            "arith.addi",
            ctx.index_t,
            _op2("arith.muli", ctx.index_t, lane, c[layout.cells_per_block]),
            cell_local,
        )
        # All threads, including padding cells, cooperate in the table copy.
        for name, array in unique.items():
            for offset in range(0, array.size, block_threads):
                index = _op2("arith.addi", ctx.index_t, linear_thread, _const_index(ctx, offset))
                guard = Operation.create(
                    "scf.if",
                    regions=2,
                    operands=[_cmpi(6, i1, index, _const_index(ctx, array.size))],
                )
                with InsertionPoint(guard.regions[0].blocks.append()):
                    indices = []
                    for axis, size in enumerate(array.shape):
                        stride = math.prod(array.shape[axis + 1 :])
                        value = _op2("arith.divui", ctx.index_t, index, _const_index(ctx, stride))
                        value = _op2("arith.remui", ctx.index_t, value, _const_index(ctx, size))
                        indices.append(value)
                    _memref_store(
                        _memref_load(ctx.global_val[name], indices, ctx.f64),
                        shared_tables[name],
                        indices,
                    )
                    Operation.create("scf.yield")
        ctx.global_val = {
            name: shared_tables[representative] for name, representative in aliases.items()
        }
        # Each coefficient is loaded once per cell, with cells contiguous in shared memory.
        for offset in range(0, layout.coefficient_size, layout.threads_per_cell):
            index = _op2("arith.addi", ctx.index_t, lane, c[offset])
            valid = _op2(
                "arith.andi", i1, valid_cell, _cmpi(6, i1, index, c[layout.coefficient_size])
            )
            guard = Operation.create("scf.if", operands=[valid], regions=2)
            with InsertionPoint(guard.regions[0].blocks.append()):
                pos = _op2(
                    "arith.addi",
                    ctx.index_t,
                    _op2("arith.muli", ctx.index_t, cell, c[layout.coefficient_size]),
                    index,
                )
                _memref_store(
                    _memref_load(coefficients, [pos], ctx.f64), shared_coeffs, [index, cell_local]
                )
                Operation.create("scf.yield")
        Operation.create("gpu.barrier")
        valid_q = _op2(
            "arith.andi", i1, valid_cell, _cmpi(6, i1, lane, c[layout.quadrature_points])
        )
        guard = Operation.create("scf.if", operands=[valid_q], regions=2)
        with InsertionPoint(guard.regions[0].blocks.append()):
            for component in range(6):
                pos = _op2(
                    "arith.addi",
                    ctx.index_t,
                    _op2("arith.muli", ctx.index_t, cell, c[6]),
                    c[component],
                )
                ctx.geometry_components[component] = _memref_load(metric, [pos], ctx.f64)
            ctx.coeffs_val = Operation.create(
                "memref.reinterpret_cast",
                results=[local_coeff],
                operands=[shared_coeffs, cell_local],
                attributes={
                    "static_offsets": DenseI64ArrayAttr.get([ShapedType.get_dynamic_size()]),
                    "static_sizes": DenseI64ArrayAttr.get([layout.coefficient_size]),
                    "static_strides": DenseI64ArrayAttr.get([layout.cells_per_block]),
                    "operandSegmentSizes": DenseI32ArrayAttr.get([1, 1, 0, 0]),
                },
            ).results[0]
            ctx.index_vars[loop_vars[0]] = lane
            cache: dict[GraphNode, Value] = {}
            for node in topo:
                if levels[node] < 2:
                    if isinstance(node, FunctionCall):
                        _, inputs, body = functions[node.function]
                        coefficient_loop = body.statements[1]
                        total = ctx.zero_f64
                        for dof_index in range(coefficient_loop.end):
                            dof = _const_index(ctx, dof_index)
                            has_point = len(inputs) > 1
                            term = _emit_coefficient_expr(
                                coefficient_loop.body.body,
                                coefficient_loop.variable,
                                dof,
                                inputs[1]._variable if has_point else None,
                                ctx.resolve_index(node.inputs[1]) if has_point else None,
                                ctx.coeffs_val,
                                table_types,
                                dict(ctx.global_val),
                                ctx.f64,
                                ctx.index_t,
                            )
                            total = _op2("arith.addf", ctx.f64, total, term)
                        cache[node] = total
                    else:
                        _emit_node(node, cache, ctx, use_signature_cache=False)
            for slot, node in enumerate(frontier):
                _memref_store(cache[node], scratch, [c[slot], lane, cell_local])
            del ctx.index_vars[loop_vars[0]]
            Operation.create("scf.yield")
        # This barrier is outside the guards: no threads exit early.
        Operation.create("gpu.barrier")
        valid_dof = _op2("arith.andi", i1, valid_cell, _cmpi(6, i1, lane, c[layout.ndofs]))
        guard = Operation.create("scf.if", operands=[valid_dof], regions=2)
        with InsertionPoint(guard.regions[0].blocks.append()):
            ctx.index_vars[loop_vars[1]] = lane
            accumulator = ctx.zero_f64
            for q_index in range(layout.quadrature_points):
                q = _const_index(ctx, q_index)
                ctx.index_vars[loop_vars[0]] = q
                cache = {
                    node: _memref_load(scratch, [c[slot], q, cell_local], ctx.f64)
                    for slot, node in enumerate(frontier)
                }
                for node in topo:
                    if levels[node] == 2:
                        _emit_node(node, cache, ctx, use_signature_cache=False)
                accumulator = _op2("arith.addf", ctx.f64, accumulator, cache[add_node.body])
                del ctx.index_vars[loop_vars[0]]
            pos = _op2(
                "arith.addi",
                ctx.index_t,
                _op2("arith.muli", ctx.index_t, cell, c[layout.ndofs]),
                lane,
            )
            row = _op1(
                "arith.index_cast",
                ctx.index_t,
                _memref_load(dofmap, [pos], IntegerType.get_signless(32)),
            )
            Operation.create(
                "memref.atomic_rmw",
                results=[ctx.f64],
                operands=[accumulator, output, row],
                attributes={"kind": IntegerAttr.get(IntegerType.get_signless(64), 0)},
            )
            del ctx.index_vars[loop_vars[1]]
            Operation.create("scf.yield")
        gpu_d.ReturnOp(operands_=[])
