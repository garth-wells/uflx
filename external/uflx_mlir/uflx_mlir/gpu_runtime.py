"""CUDA/HIP device-memory execution for generated linear-form assembly kernels."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from mlir.ir import Module

from uflx_mlir.geometry import geometry_kernel_name
from uflx_mlir.gpu_assembly import (
    assemble_amdgcn_to_hsaco,
    extract_amdgcn_text,
    extract_ptx_text,
    lower_module_to_nvvm,
    lower_module_to_rocdl,
)
from uflx_mlir.gpu_driver import _launch, launch_isolated
from uflx_mlir.gpu_linear import LinearAssemblyLayout


def assemble_linear_gpu(
    module: Module,
    layout: LinearAssemblyLayout,
    kernel_name: str,
    coordinates: np.ndarray,
    coefficients: np.ndarray,
    cell_dofs: np.ndarray,
    output: np.ndarray,
    *,
    backend: str = "cuda",
    chip: str = "sm_80",
    device: int = 0,
    rocm_path: str = "/opt/rocm",
    geometry: np.ndarray | None = None,
) -> float:
    """Accumulate into output using CUDA or HIP device allocations.

    Inputs use the shapes documented by generate_linear_assembly_gpu_module.
    The module is lowered in place; generate a fresh module for each call.
    Nonzero initial output is preserved. Returned seconds measure one launch
    and synchronization, excluding compilation, allocations and transfers.
    Stored geometry is computed by a separate GPU kernel before timing unless
    supplied explicitly. No explicit action warm-up is performed. Empty meshes
    return zero without a launch.
    HIP execution uses a fresh worker process to isolate ROCm LLVM from MLIR.
    Worker startup and array exchange are excluded from the returned timing.

    Args:
        module: Fresh generated GPU module, consumed by lowering.
        layout: Generator-provided input and launch layout.
        kernel_name: Generated kernel symbol.
        coordinates: Contiguous f64 array (ncells, coordinate DOFs, dimension).
        coefficients: Contiguous f64 array (ncells, packed coefficient size).
        cell_dofs: Contiguous i32 array (ncells, test DOFs).
        output: Writable contiguous f64 global vector, updated in place.
        backend: "cuda" or "amd" (HIP).
        chip: Architecture, for example "sm_89" or "gfx1100".
        device: Device ordinal within the selected backend.
        rocm_path: ROCm installation used for AMD compilation and runtime.
        geometry: Optional C-contiguous float64 array (ncells, layout.geometry_size)
            containing abs(det(J)) * inv(J) * inv(J).T in packed symmetric order,
            without quadrature weights. Supplying it skips geometry preparation.
            Only supported when layout.geometry_size is nonzero.

    Returns:
        Synchronized kernel-launch wall time in seconds.

    Raises:
        ValueError: Invalid arrays or launch configuration.
        RuntimeError: A driver operation fails.
    """
    if backend not in ("cuda", "amd"):
        raise ValueError("backend must be 'cuda' or 'amd'")
    if coordinates.ndim != 3 or coordinates.shape[1:] != layout.coordinate_shape:
        raise ValueError("coordinates must have shape (ncells, *layout.coordinate_shape)")
    ncells = coordinates.shape[0]
    if coefficients.shape != (ncells, layout.coefficient_size):
        raise ValueError("coefficients have the wrong packed per-cell shape")
    if cell_dofs.shape != (ncells, layout.ndofs):
        raise ValueError("cell_dofs have the wrong per-cell shape")
    if output.ndim != 1 or not output.flags.writeable:
        raise ValueError("output must be a writable vector")
    for array, dtype in (
        (output, np.float64),
        (coordinates, np.float64),
        (coefficients, np.float64),
        (cell_dofs, np.int32),
    ):
        if array.dtype != dtype or not array.flags.c_contiguous:
            raise ValueError(f"arrays must be C-contiguous with dtype {dtype}")
    if geometry is not None:
        if not layout.geometry_size:
            raise ValueError("this kernel does not accept stored geometry")
        if geometry.shape != (ncells, layout.geometry_size):
            raise ValueError("geometry has the wrong per-cell shape")
        if geometry.dtype != np.float64 or not geometry.flags.c_contiguous:
            raise ValueError("geometry must be C-contiguous with dtype float64")
    if ncells == 0:
        return 0.0
    if cell_dofs.min() < 0 or cell_dofs.max() >= len(output):
        raise ValueError("cell_dofs contain an out-of-range global DOF")
    geometry_coordinates = None
    geometry_input = coordinates
    if layout.geometry_size:
        if geometry is None:
            geometry = np.zeros((ncells, layout.geometry_size), dtype=np.float64)
            geometry_coordinates = coordinates
        geometry_input = geometry
    rocm = Path(rocm_path)
    if backend == "cuda":
        lower_module_to_nvvm(module, cubin_chip=chip)
        code = extract_ptx_text(module).encode()
        library = "libcuda.so.1"
    else:
        lower_module_to_rocdl(module, chip=chip, link_device_libraries=False)
        code = assemble_amdgcn_to_hsaco(
            extract_amdgcn_text(module), chip=chip, toolkit_path=str(rocm)
        )
        library = str(rocm / "lib/libamdhip64.so")
    launch = launch_isolated if backend == "amd" else _launch
    return launch(
        library,
        code,
        kernel_name,
        [output, geometry_input.reshape(-1), coefficients.reshape(-1), cell_dofs.reshape(-1)],
        layout.grid_shape(ncells),
        layout.block_shape,
        backend,
        device,
        geometry_coordinates=geometry_coordinates,
        geometry_name=geometry_kernel_name(kernel_name)
        if geometry_coordinates is not None
        else None,
    )
