"""GPU driver execution and HIP worker; deliberately independent of MLIR imports.

Run this file by path in a fresh interpreter, not via the uflx_mlir package.
"""

from __future__ import annotations

import ctypes as ct
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np


def _launch(
    library,
    code,
    name,
    arrays,
    grid,
    block,
    backend,
    device,
    *,
    geometry_coordinates=None,
    geometry_name=None,
) -> float:
    """Launch rank-1 memrefs; release all allocated resources on errors."""
    if geometry_coordinates is not None and geometry_name is None:
        raise ValueError("geometry_name is required for geometry setup")
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
                assert geometry_name is not None
                setup = pointer()
                check(
                    get_function(ct.byref(setup), module, geometry_name.encode()),
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


def _worker(directory: Path) -> None:
    """Execute a serialized HIP request without importing the compiler package."""
    request = json.loads((directory / "request.json").read_text())
    arrays = [
        np.load(directory / f"array{i}.npy", mmap_mode="r+" if i == 0 else "r", allow_pickle=False)
        for i in range(request["array_count"])
    ]
    coordinates = (
        np.load(directory / "coordinates.npy", mmap_mode="r", allow_pickle=False)
        if request["geometry_name"] is not None
        else None
    )
    elapsed = _launch(
        request["library"],
        (directory / "kernel.hsaco").read_bytes(),
        request["name"],
        arrays,
        request["grid"],
        request["block"],
        "amd",
        request["device"],
        geometry_coordinates=coordinates,
        geometry_name=request["geometry_name"],
    )
    arrays[0].flush()
    (directory / "result.json").write_text(json.dumps({"elapsed": elapsed}))


def launch_isolated(
    library,
    code,
    name,
    arrays,
    grid,
    block,
    backend,
    device,
    *,
    geometry_coordinates=None,
    geometry_name=None,
) -> float:
    """Run HIP in a fresh interpreter, keeping ROCm LLVM separate from MLIR LLVM."""
    if backend != "amd":
        raise ValueError("The isolated worker supports only HIP")
    with tempfile.TemporaryDirectory(prefix="uflx-hip-worker-") as temporary:
        directory = Path(temporary)
        (directory / "kernel.hsaco").write_bytes(code)
        for i, array in enumerate(arrays):
            np.save(directory / f"array{i}.npy", array, allow_pickle=False)
        if geometry_coordinates is not None:
            np.save(directory / "coordinates.npy", geometry_coordinates, allow_pickle=False)
        request = dict(
            array_count=len(arrays),
            library=library,
            name=name,
            grid=grid,
            block=block,
            device=device,
            geometry_name=geometry_name,
        )
        (directory / "request.json").write_text(json.dumps(request))
        worker = Path(__file__).with_name("gpu_driver.py")
        # Execute by path: importing this through the package would import MLIR.
        result = subprocess.run(
            [sys.executable, "-I", os.fspath(worker), os.fspath(directory)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode:
            raise RuntimeError(
                f"HIP worker exited with status {result.returncode}:\n"
                f"{result.stderr or result.stdout}"
            )
        elapsed = float(json.loads((directory / "result.json").read_text())["elapsed"])
        # Publish output only after a successful launch and worker shutdown.
        np.copyto(arrays[0], np.load(directory / "array0.npy", allow_pickle=False))
        return elapsed


if __name__ == "__main__":
    _worker(Path(sys.argv[1]))
