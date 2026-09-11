"""MLIR code-generation backend for UFLx forms, built on the Python op-builder API.

See README.md for what this package does, its current restrictions, and
how to build the `mlir` Python bindings it depends on (not a pip
dependency -- see below).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from uflx_mlir.geometry import geometry_kernel_name

if TYPE_CHECKING:
    from uflx_mlir.emit import generate_mlir_module

__all__ = ["generate_mlir_module", "geometry_kernel_name"]


def __getattr__(name: str) -> Any:
    """Load the MLIR-dependent generator only when it is requested."""
    if name == "generate_mlir_module":
        from uflx_mlir.emit import generate_mlir_module

        return generate_mlir_module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
