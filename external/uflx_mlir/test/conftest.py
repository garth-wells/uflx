"""Shared capability checks for the MLIR backend tests."""

import pytest


def _require_gpu_target(target: str) -> None:
    """Probe the loaded bindings, skipping only an explicitly missing LLVM target."""
    from mlir.ir import Context, Location, MLIRError, Module

    from uflx_mlir.gpu_assembly import lower_module_to_nvvm, lower_module_to_rocdl

    with Context(), Location.unknown():
        module = Module.parse(
            "module attributes {gpu.container_module} {"
            "  gpu.module @target_probe {"
            "    gpu.func @probe() kernel { gpu.return }"
            "  }"
            "}"
        )
        try:
            if target == "AMDGPU":
                lower_module_to_rocdl(module, chip="gfx90a", link_device_libraries=False)
            else:
                assert target == "NVPTX"
                lower_module_to_nvvm(module, cubin_chip="sm_80")
        except MLIRError as error:
            if f"The `{target}` target was not built" in str(error):
                pytest.skip(f"MLIR was built without the {target} target")
            raise


@pytest.fixture(scope="session")
def require_amdgpu() -> None:
    """Require AMD output support independently of ROCm tools or GPU availability."""
    _require_gpu_target("AMDGPU")


@pytest.fixture(scope="session")
def require_nvptx() -> None:
    """Require NVIDIA output support independently of CUDA tools or GPU availability."""
    _require_gpu_target("NVPTX")
