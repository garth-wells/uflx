"""Batched GPU vector assembly for real linear forms with packed coefficients."""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast

import basix
import numpy as np
from mlir.dialects import gpu as gpu_d
from mlir.ir import (
    Attribute,
    Context,
    DenseElementsAttr,
    DenseI32ArrayAttr,
    DenseI64ArrayAttr,
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
    StridedLayoutAttr,
    StringAttr,
    Type,
    TypeAttr,
    UnitAttr,
    Value,
)
from mlir.passmanager import PassManager
from uflx.geometry import CoordinateDofComponent
from uflx.graphs import as_graph
from uflx_codegeneration import symbols
from uflx_codegeneration.nodes import ArrayEntry, Block, Loop
from uflx_codegeneration.quadrature import QuadratureLoop

from uflx_mlir.emit import (
    _build_nest,
    _const_f64,
    _const_index,
    _emit_affine_tetrahedron_geometry_values,
    _emit_coefficient_function,
    _emit_fission_scratch_allocas,
    _memref_load,
    _memref_store,
    _op1,
    _op2,
    _OpCtx,
)
from uflx_mlir.geometry import geometry_kernel_name
from uflx_mlir.gpu_assembly import _cmpi, gpu_module_name
from uflx_mlir.gpu_cooperative import cooperative_memory_bytes, emit_cooperative_action
from uflx_mlir.hoist import (
    compute_fission_plan,
    distribute_shallow_factors,
    reorder_quadrature_outermost,
    topo_order,
    walk_loop_chain,
)
from uflx_mlir.lowering import collect_int_constants, coordinate_shape, lower_form


@dataclass(frozen=True)
class LinearAssemblyLayout:
    """Input strides and launch dimensions for one generated vector kernel."""

    ndofs: int
    quadrature_points: int
    coordinate_shape: tuple[int, int]
    coefficient_size: int
    threads_per_cell: int
    cells_per_block: int
    geometry_size: int = 0
    cooperative: bool = False

    @property
    def block_shape(self) -> tuple[int, int, int]:
        """Return cell/point dimensions for cooperative or serial quadrature."""
        if self.cooperative:
            return self.cells_per_block, self.threads_per_cell, 1
        return self.threads_per_cell, 1, self.cells_per_block

    def grid_shape(self, ncells: int) -> tuple[int, int, int]:
        """Round up the grid to include a partially populated final block."""
        if ncells < 0:
            raise ValueError("ncells must be nonnegative")
        return (ncells + self.cells_per_block - 1) // self.cells_per_block, 1, 1


def _coefficient_size(functions) -> int:
    """Recover the CPU lowering's packed per-cell coefficient stride."""
    size = 0
    for _, _, body in functions.values():
        assert isinstance(body, Block)
        loop = body.statements[1]
        assert isinstance(loop, Loop) and isinstance(loop.end, int)
        for node in as_graph(loop.body):
            if isinstance(node, ArrayEntry) and node.array == symbols.coefficients:
                index = str(node.index[0])
                match = re.fullmatch(r"(?:(\d+)\s*\+\s*)?" + re.escape(loop.variable), index)
                if match is None:
                    raise NotImplementedError(f"Unsupported coefficient offset: {index}")
                size = max(size, int(match[1] or 0) + loop.end)
    return size


def _lower_device_abs(operation) -> None:
    """Lower f64 absolute value without LLVM fabs intrinsics unsupported by NVVM."""
    for region in operation.regions:
        for block in region.blocks:
            for child in list(block.operations):
                _lower_device_abs(child.operation)
                if child.operation.name != "math.absf":
                    continue
                with InsertionPoint(child):
                    i64 = IntegerType.get_signless(64)
                    bits = _op1("arith.bitcast", i64, child.operands[0])
                    mask = Operation.create(
                        "arith.constant",
                        results=[i64],
                        attributes={"value": IntegerAttr.get(i64, (1 << 63) - 1)},
                    ).results[0]
                    magnitude = _op2("arith.andi", i64, bits, mask)
                    result = _op1("arith.bitcast", child.results[0].type, magnitude)
                    child.results[0].replace_all_uses_with(result)
                child.operation.erase()


def _gpu_kernel(name: str, arg_types: list[Type]) -> Operation:
    """Create a kernel without relying on version-specific generated constructors."""
    return Operation.create(
        "gpu.func",
        regions=1,
        attributes={
            "sym_name": StringAttr.get(name),
            "function_type": TypeAttr.get(FunctionType.get(arg_types, [])),
            "gpu.kernel": UnitAttr.get(),
        },
    )


def _emit_geometry_gpu_kernel(kernel_name, flat, f64, index_t, i1, coords_shape):
    """Store six affine metric components per cell in blocks of 128 x threads."""
    kernel = _gpu_kernel(geometry_kernel_name(kernel_name), [flat, flat])
    entry = kernel.regions[0].blocks.append(flat, flat)
    ctx = _OpCtx((6,), coords_shape, {}, f64, index_t)
    with InsertionPoint(entry):
        metric, coords = cast(Sequence[Value], entry.arguments)
        c = {v: _const_index(ctx, v) for v in range(13)}
        cell_id = _op2(
            "arith.addi",
            index_t,
            _op2(
                "arith.muli",
                index_t,
                gpu_d.block_id(gpu_d.Dimension.x),
                _const_index(ctx, 128),
            ),
            gpu_d.thread_id(gpu_d.Dimension.x),
        )
        length = Operation.create("memref.dim", results=[index_t], operands=[coords, c[0]]).results[
            0
        ]
        ncells = _op2("arith.divui", index_t, length, c[12])
        guard = Operation.create("scf.if", operands=[_cmpi(6, i1, cell_id, ncells)], regions=2)
        with InsertionPoint(guard.regions[0].blocks.append()):
            local = Operation.create(
                "memref.alloca", results=[MemRefType.get(list(coords_shape), f64)]
            ).results[0]
            base = _op2("arith.muli", index_t, cell_id, c[12])
            for point in range(4):
                for component in range(3):
                    pos = _op2("arith.addi", index_t, base, c[3 * point + component])
                    _memref_store(_memref_load(coords, [pos], f64), local, [c[point], c[component]])
            values = _emit_affine_tetrahedron_geometry_values(local, c, f64)
            base = _op2("arith.muli", index_t, cell_id, c[6])
            for component, value in enumerate(values):
                pos = _op2("arith.addi", index_t, base, c[component])
                _memref_store(value, metric, [pos])
            Operation.create("scf.yield")
        gpu_d.ReturnOp(operands_=[])


def generate_linear_assembly_gpu_module(
    form,
    degree: int,
    kernel_name: str,
    cell: basix.CellType,
    *,
    cells_per_block: int | None = None,
    target_block_size: int | None = None,
    precompute_geometry: bool = True,
    cooperative: bool = True,
) -> tuple[Module, LinearAssemblyLayout]:
    """Generate a GPU kernel that accumulates a linear form into a global vector.

    The device ABI is four flattened rank-1 memrefs: output f64 vector,
    per-cell geometry (or coordinate) f64 values, packed coefficient f64 values,
    and per-cell test-space i32 global DOF indices. Coordinates have shape
    (ncells, *layout.coordinate_shape); coefficients have shape
    (ncells, layout.coefficient_size). Coefficients within each cell follow
    the CPU emitter's order (increasing coefficient count). With no
    coefficients, pass an empty buffer. All buffers must reside on device.
    Output is additive: zero it first for fresh assembly.

    Affine stiffness forms use a stored metric by default: the second buffer
    has shape (ncells, 6), in order (G00, G01, G02, G11, G12, G22), with
    G = abs(det(J)) * inv(J) * inv(J).T. Quadrature weights remain in the action.
    A companion <kernel_name>_geometry GPU kernel fills this buffer from
    coordinates once per mesh update. Its ABI is (metric, coordinates), both
    flattened f64 memrefs; launch blocks of (128, 1, 1), one x thread per cell. Other forms retain
    the coordinate ABI. layout.geometry_size identifies the selected ABI.

    By default, supported stored-metric forms use cooperative quadrature.
    The x block dimension selects a cell, and y selects a quadrature point
    during evaluation or a test DOF during contraction. Coefficients and basis
    tables are cached in aligned shared memory. Quadrature-only values are
    evaluated once, shared across test DOFs, and contracted with one atomic
    add per output DOF. Fixed-size contractions are unrolled. All block threads
    participate in both barriers, including inactive cells in the final block.

    The serial fallback assigns one x thread per test DOF and packs cells in z.
    It is selected by cooperative=False, unsupported geometry, or shared storage
    exceeding 48 KiB. Both paths retain the portable 256-thread block limit.

    Args:
        form: Real linear form with one integral and one test-space DOF axis.
        degree: Polynomial degree used by the existing quadrature selector.
        kernel_name: Exported GPU kernel symbol.
        cell: Basix reference cell supported by the CPU lowering.
        cells_per_block: Explicit cell grouping, or None for automatic grouping.
        target_block_size: Desired automatic block size, from 1 to 256. None
            selects 256 for cooperative quadrature and 128 for serial quadrature.
        precompute_geometry: Extract supported affine metrics into a separate GPU
            setup kernel. False retains the previous inline-geometry implementation.
        cooperative: Share coefficient evaluation and quadrature fluxes across
            test DOFs when supported. False retains serial quadrature per test DOF.

    Returns:
        A module accepted by lower_module_to_nvvm/lower_module_to_rocdl,
        and a layout defining input strides and launch dimensions.

    Raises:
        ValueError: Invalid launch configuration.
        NotImplementedError: Unsupported form shape or more than 256 test DOFs.
    """
    requested_cells, requested_target = cells_per_block, target_block_size
    if target_block_size is not None and (
        not isinstance(target_block_size, int) or not 1 <= target_block_size <= 256
    ):
        raise ValueError("target_block_size must be an integer from 1 to 256")
    tables, graph, geometry, functions = lower_form(form, degree, cell)
    chain, add_node = walk_loop_chain(graph.root)
    chain = reorder_quadrature_outermost(chain)
    if (
        len(chain) != 2
        or not isinstance(chain[0][0], QuadratureLoop)
        or not isinstance(chain[1][0], Loop)
        or len(add_node.shape) != 1
    ):
        raise NotImplementedError(
            "GPU vector assembly requires one quadrature and one test DOF axis"
        )
    quadrature, _ = chain[0]
    dof_loop, dof_var = chain[1]
    assert isinstance(quadrature, QuadratureLoop) and isinstance(dof_loop, Loop)
    ndofs = dof_loop.end
    if not isinstance(ndofs, int) or dof_loop.start != 0 or not 1 <= ndofs <= 256:
        raise NotImplementedError("GPU vector assembly requires 1 to 256 test DOFs")
    coords_shape = coordinate_shape(form)
    ncoords = int(np.prod(coords_shape))
    ncoeff = _coefficient_size(functions)
    # A partially extracted form can still need coordinates (e.g. mass + stiffness).
    # Preserve its existing coordinate ABI until all geometric factors are extracted.
    stored_geometry = (
        precompute_geometry
        and geometry is not None
        and not any(isinstance(node, CoordinateDofComponent) for node in graph)
    )
    cooperative = bool(
        cooperative and stored_geometry and functions and quadrature.rule.npoints <= 256
    )
    if target_block_size is None:
        target_block_size = 256 if cooperative else 128
    threads = max(ndofs, quadrature.rule.npoints) if cooperative else 1 << (ndofs - 1).bit_length()
    if cells_per_block is None:
        groups = max(1, min(64, target_block_size // threads))
        cells_per_block = groups if cooperative else 1 << (groups.bit_length() - 1)
    if (
        not isinstance(cells_per_block, int)
        or not 1 <= cells_per_block <= 64
        or threads * cells_per_block > 256
    ):
        raise ValueError("cells_per_block must be 1..64 and total threads must not exceed 256")
    geometry_size = geometry.output_size if stored_geometry and geometry is not None else 0
    input_stride = geometry_size or ncoords
    layout = LinearAssemblyLayout(
        ndofs,
        quadrature.rule.npoints,
        coords_shape,
        ncoeff,
        threads,
        cells_per_block,
        geometry_size,
        cooperative,
    )
    loop_vars = [var for _, var in chain]
    distribute_shallow_factors(add_node, loop_vars)
    if (
        cooperative
        and cooperative_memory_bytes(add_node, chain, tables, ncoeff, cells_per_block) > 48 * 1024
    ):
        return generate_linear_assembly_gpu_module(
            form,
            degree,
            kernel_name,
            cell,
            cells_per_block=requested_cells,
            target_block_size=requested_target,
            precompute_geometry=precompute_geometry,
            cooperative=False,
        )
    constants = collect_int_constants(graph.root, add_node.shape)
    constants.update(
        [ncoords, ncoeff, input_stride, cells_per_block, *range(max(coords_shape) + 1)]
    )
    if geometry is not None:
        constants.update(range(geometry.output_size))

    with Context(), Location.unknown():
        module = Module.create()
        module.operation.attributes["gpu.container_module"] = UnitAttr.get()
        f64, index_t, i32, i1 = (
            F64Type.get(),
            IndexType.get(),
            IntegerType.get_signless(32),
            IntegerType.get_signless(1),
        )
        dyn = ShapedType.get_dynamic_size()
        flat = MemRefType.get([dyn], f64)
        indices = MemRefType.get([dyn], i32)
        local_coeff = MemRefType.get([ncoeff], f64, StridedLayoutAttr.get(dyn, [1]))
        if cooperative:
            local_coeff = MemRefType.get(
                [ncoeff],
                f64,
                StridedLayoutAttr.get(dyn, [cells_per_block]),
                memory_space=Attribute.parse("#gpu.address_space<workgroup>"),
            )
        ctx = _OpCtx(
            add_node.shape, coords_shape, {k: v.shape for k, v in tables.items()}, f64, index_t
        )
        with InsertionPoint(module.body):
            gmod = Operation.create(
                "gpu.module",
                regions=1,
                attributes={"sym_name": StringAttr.get(gpu_module_name(kernel_name))},
            )
            body = gmod.regions[0].blocks.append()
            with InsertionPoint(body):
                table_types = {}
                for name, arr in sorted(tables.items()):
                    table_types[name] = MemRefType.get(list(arr.shape), f64)
                    Operation.create(
                        "memref.global",
                        attributes={
                            "sym_name": StringAttr.get(name),
                            "sym_visibility": StringAttr.get("private"),
                            "type": TypeAttr.get(table_types[name]),
                            "constant": UnitAttr.get(),
                            "initial_value": DenseElementsAttr.get(
                                np.ascontiguousarray(arr),
                                type=RankedTensorType.get(list(arr.shape), f64),
                            ),
                        },
                    )
                for name, (_, inputs, function_body) in (
                    sorted(functions.items()) if not cooperative else []
                ):
                    _emit_coefficient_function(
                        name, inputs, function_body, table_types, f64, index_t, local_coeff
                    )
                for op in body.operations:
                    if op.operation.name == "func.func":
                        op.attributes["sym_visibility"] = StringAttr.get("private")
                if stored_geometry:
                    _emit_geometry_gpu_kernel(kernel_name, flat, f64, index_t, i1, coords_shape)
                arg_types: list[Type] = [flat, flat, flat, indices]
                kernel = _gpu_kernel(kernel_name, arg_types)
                if cooperative:
                    emit_cooperative_action(
                        kernel,
                        arg_types,
                        ctx,
                        layout,
                        chain,
                        add_node,
                        table_types,
                        local_coeff,
                        tables,
                        functions,
                    )
                else:
                    entry = kernel.regions[0].blocks.append(*arg_types)
                    with InsertionPoint(entry):
                        output, coords, coefficients, dofmap = cast(
                            Sequence[Value], entry.arguments
                        )
                        for value in sorted(constants):
                            ctx.index_const[value] = _const_index(ctx, value)
                        ctx.zero_f64 = _const_f64(ctx, 0.0)
                        c = ctx.index_const
                        x = cast(Value, gpu_d.thread_id(gpu_d.Dimension.x))
                        z = gpu_d.thread_id(gpu_d.Dimension.z)
                        cell_id = _op2(
                            "arith.addi",
                            index_t,
                            _op2(
                                "arith.muli",
                                index_t,
                                gpu_d.block_id(gpu_d.Dimension.x),
                                c[cells_per_block],
                            ),
                            z,
                        )
                        length = Operation.create(
                            "memref.dim", results=[index_t], operands=[coords, c[0]]
                        ).results[0]
                        ncells = _op2("arith.divui", index_t, length, c[input_stride])
                        valid = _op2(
                            "arith.andi",
                            i1,
                            _cmpi(6, i1, cell_id, ncells),
                            _cmpi(6, i1, x, c[ndofs]),
                        )
                        guard = Operation.create("scf.if", operands=[valid], regions=2)
                        with InsertionPoint(guard.regions[0].blocks.append()):
                            for name, ty in table_types.items():
                                ctx.global_val[name] = Operation.create(
                                    "memref.get_global",
                                    results=[ty],
                                    attributes={"name": FlatSymbolRefAttr.get(name)},
                                ).results[0]
                            if stored_geometry:
                                base = _op2("arith.muli", index_t, cell_id, c[geometry_size])
                                for component in range(geometry_size):
                                    pos = _op2("arith.addi", index_t, base, c[component])
                                    ctx.geometry_components[component] = _memref_load(
                                        coords, [pos], f64
                                    )
                            else:
                                local_coords = Operation.create(
                                    "memref.alloca",
                                    results=[MemRefType.get(list(coords_shape), f64)],
                                ).results[0]
                                base = _op2("arith.muli", index_t, cell_id, c[ncoords])
                                for point in range(coords_shape[0]):
                                    for component in range(coords_shape[1]):
                                        offset = _const_index(
                                            ctx, point * coords_shape[1] + component
                                        )
                                        pos = _op2("arith.addi", index_t, base, offset)
                                        _memref_store(
                                            _memref_load(coords, [pos], f64),
                                            local_coords,
                                            [c[point], c[component]],
                                        )
                                ctx.coords_val = local_coords  # type: ignore[attr-defined]
                                if geometry is not None:
                                    ctx.geometry_components.update(
                                        enumerate(
                                            _emit_affine_tetrahedron_geometry_values(
                                                local_coords, c, f64
                                            )
                                        )
                                    )
                            if functions:
                                offset = _op2("arith.muli", index_t, cell_id, c[ncoeff])
                                ctx.coeffs_val = Operation.create(
                                    "memref.reinterpret_cast",
                                    results=[local_coeff],
                                    operands=[coefficients, offset],
                                    attributes={
                                        "static_offsets": DenseI64ArrayAttr.get([dyn]),
                                        "static_sizes": DenseI64ArrayAttr.get([ncoeff]),
                                        "static_strides": DenseI64ArrayAttr.get([1]),
                                        "operandSegmentSizes": DenseI32ArrayAttr.get([1, 1, 0, 0]),
                                    },
                                ).results[0]
                            ctx.thread_bindings = {dof_var: x}
                            accum = Operation.create(
                                "memref.alloca", results=[MemRefType.get([], f64)]
                            ).results[0]
                            _memref_store(ctx.zero_f64, accum, [])

                            def commit(context, result):
                                _memref_store(
                                    _op2("arith.addf", f64, _memref_load(accum, [], f64), result),
                                    accum,
                                    [],
                                )

                            ctx.commit = commit
                            levels, groups = compute_fission_plan(add_node, loop_vars)
                            _emit_fission_scratch_allocas(groups, chain, ctx)
                            by_depth = {}
                            for group in groups:
                                by_depth.setdefault(group.depth, []).append(group)
                            _build_nest(
                                0,
                                chain,
                                levels,
                                topo_order(add_node),
                                set(),
                                {},
                                ctx,
                                add_node,
                                by_depth,
                            )
                            pos = _op2(
                                "arith.addi",
                                index_t,
                                _op2("arith.muli", index_t, cell_id, c[ndofs]),
                                x,
                            )
                            row = _op1(
                                "arith.index_cast", index_t, _memref_load(dofmap, [pos], i32)
                            )
                            Operation.create(
                                "memref.atomic_rmw",
                                results=[f64],
                                operands=[_memref_load(accum, [], f64), output, row],
                                attributes={
                                    "kind": IntegerAttr.get(IntegerType.get_signless(64), 0)
                                },
                            )
                            Operation.create("scf.yield")
                        gpu_d.ReturnOp(operands_=[])
                try:
                    gpu_d.ModuleEndOp()
                except AttributeError:
                    pass
        module.operation.verify()
        # Inline the CPU lowering's coefficient helpers into the GPU kernel so
        # no host function calls or host table addresses reach the device code.
        PassManager.parse("builtin.module(gpu.module(inline,canonicalize,cse))").run(
            module.operation
        )
        _lower_device_abs(module.operation)
        module.operation.verify()
        return module, layout
