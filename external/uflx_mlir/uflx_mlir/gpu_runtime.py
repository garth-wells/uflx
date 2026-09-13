"""CUDA/HIP device-memory execution for generated linear-form assembly kernels."""

from __future__ import annotations

import ctypes as ct
import os
import tempfile
import time
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
    return _launch(
        library,
        code,
        kernel_name,
        [output, geometry_input.reshape(-1), coefficients.reshape(-1), cell_dofs.reshape(-1)],
        layout.grid_shape(ncells),
        layout.block_shape,
        backend,
        device,
        geometry_coordinates=geometry_coordinates,
    )


def _launch(
    library, code, name, arrays, grid, block, backend, device, *, geometry_coordinates=None
) -> float:
    """Launch four rank-1 memrefs; release all allocated resources on errors."""
    driver = ct.CDLL(library)
    pointer = ct.c_void_p
    cuda = backend == "cuda"
    address_type = ct.c_uint64 if cuda else pointer

    def bind(name, *args):
        fn = getattr(driver, name)
        fn.restype = ct.c_int
        fn.argtypes = list(args)
        return fn

    init = bind("cuInit" if cuda else "hipInit", ct.c_uint)
    get_function = bind(
        "cuModuleGetFunction" if cuda else "hipModuleGetFunction",
        ct.POINTER(pointer),
        pointer,
        ct.c_char_p,
    )
    allocate = bind("cuMemAlloc_v2" if cuda else "hipMalloc", ct.POINTER(address_type), ct.c_size_t)
    free = bind("cuMemFree_v2" if cuda else "hipFree", address_type)
    unload = bind("cuModuleUnload" if cuda else "hipModuleUnload", pointer)
    launch = bind(
        "cuLaunchKernel" if cuda else "hipModuleLaunchKernel",
        pointer,
        ct.c_uint,
        ct.c_uint,
        ct.c_uint,
        ct.c_uint,
        ct.c_uint,
        ct.c_uint,
        ct.c_uint,
        pointer,
        ct.POINTER(pointer),
        ct.POINTER(pointer),
    )
    sync = bind("cuCtxSynchronize" if cuda else "hipDeviceSynchronize")
    api = {}
    if cuda:
        get_error = bind("cuGetErrorString", ct.c_int, ct.POINTER(ct.c_char_p))
        load = bind("cuModuleLoadData", ct.POINTER(pointer), pointer)
        api["htod"] = bind("cuMemcpyHtoD_v2", address_type, pointer, ct.c_size_t)
        api["dtoh"] = bind("cuMemcpyDtoH_v2", pointer, address_type, ct.c_size_t)
        get_device = bind("cuDeviceGet", ct.POINTER(ct.c_int), ct.c_int)
        api["retain"] = bind("cuDevicePrimaryCtxRetain", ct.POINTER(pointer), ct.c_int)
        api["release"] = bind("cuDevicePrimaryCtxRelease_v2", ct.c_int)
        api["push"] = bind("cuCtxPushCurrent_v2", pointer)
        api["pop"] = bind("cuCtxPopCurrent_v2", ct.POINTER(pointer))
    else:
        get_error = driver.hipGetErrorString
        get_error.argtypes = [ct.c_int]
        get_error.restype = ct.c_char_p
        load = bind("hipModuleLoad", ct.POINTER(pointer), ct.c_char_p)
        api["copy"] = bind("hipMemcpy", pointer, pointer, ct.c_size_t, ct.c_int)
        api["set_device"] = bind("hipSetDevice", ct.c_int)
        get_device = bind("hipGetDevice", ct.POINTER(ct.c_int))

    def check(status, operation):
        if status:
            if cuda:
                message = ct.c_char_p()
                get_error(status, ct.byref(message))
                detail = message.value
            else:
                detail = get_error(status)
            raise RuntimeError(f"{operation}: {detail.decode() if detail else status}")

    check(init(0), "initialize GPU driver")
    module, function, context = pointer(), pointer(), pointer()
    allocations = []
    previous = ct.c_int()
    selected = ct.c_int()
    retained = pushed = selected_hip = False
    cleanup = []
    try:
        if cuda:
            check(get_device(ct.byref(selected), device), "cuDeviceGet")
            check(api["retain"](ct.byref(context), selected), "cuDevicePrimaryCtxRetain")
            retained = True
            check(api["push"](context), "cuCtxPushCurrent")
            pushed = True
        else:
            check(get_device(ct.byref(previous)), "hipGetDevice")
            check(api["set_device"](device), "hipSetDevice")
            selected_hip = True
        with tempfile.TemporaryDirectory(prefix="uflx-linear-") as directory:
            if cuda:
                blob = ct.create_string_buffer(code)
                check(load(ct.byref(module), blob), "cuModuleLoadData")
            else:
                path = Path(directory) / "kernel.hsaco"
                path.write_bytes(code)
                check(load(ct.byref(module), os.fsencode(path)), "hipModuleLoad")
            check(get_function(ct.byref(function), module, name.encode()), "get kernel function")

            def upload(array):
                """Upload an array and return its rank-1 memref arguments."""
                address = address_type()
                check(allocate(ct.byref(address), max(1, array.nbytes)), "allocate device array")
                allocations.append(address)
                if array.nbytes:
                    status = (
                        api["htod"](address, pointer(array.ctypes.data), array.nbytes)
                        if cuda
                        else api["copy"](address, pointer(array.ctypes.data), array.nbytes, 1)
                    )
                    check(status, "copy to device")
                address_value = address.value
                if address_value is None:
                    raise RuntimeError("Device allocation returned a null pointer")
                return [
                    address_type(address_value),
                    address_type(address_value),
                    ct.c_int64(0),
                    ct.c_int64(array.size),
                    ct.c_int64(1),
                ]

            arguments = []
            for array in arrays:
                arguments.extend(upload(array))
            if geometry_coordinates is not None:
                setup = pointer()
                check(
                    get_function(ct.byref(setup), module, geometry_kernel_name(name).encode()),
                    "get geometry kernel function",
                )
                geometry_args = arguments[5:10] + upload(geometry_coordinates.reshape(-1))
                geometry_params = (pointer * len(geometry_args))(
                    *(ct.cast(ct.byref(x), pointer) for x in geometry_args)
                )
                ncells = geometry_coordinates.shape[0]
                check(
                    launch(
                        setup,
                        (ncells + 127) // 128,
                        1,
                        1,
                        128,
                        1,
                        1,
                        0,
                        None,
                        geometry_params,
                        None,
                    ),
                    "precompute geometry",
                )
                check(sync(), "synchronize geometry setup")
            parameters = (pointer * len(arguments))(
                *(ct.cast(ct.byref(x), pointer) for x in arguments)
            )
            start = time.perf_counter()
            check(launch(function, *grid, *block, 0, None, parameters, None), "launch kernel")
            check(sync(), "synchronize kernel")
            elapsed = time.perf_counter() - start
            output = arrays[0]
            status = (
                api["dtoh"](pointer(output.ctypes.data), allocations[0], output.nbytes)
                if cuda
                else api["copy"](pointer(output.ctypes.data), allocations[0], output.nbytes, 2)
            )
            check(status, "copy output to host")
    finally:
        import sys

        failed = sys.exc_info()[0] is not None
        for address in allocations:
            cleanup.append((free(address), "free device array"))
        if module.value:
            cleanup.append((unload(module), "unload module"))
        if cuda:
            if pushed:
                restored = pointer()
                cleanup.append((api["pop"](ct.byref(restored)), "restore CUDA context"))
            if retained:
                cleanup.append((api["release"](selected), "release CUDA context"))
        elif selected_hip:
            cleanup.append((api["set_device"](previous), "restore HIP device"))
        if not failed:
            for status, operation in cleanup:
                check(status, operation)
    return elapsed
