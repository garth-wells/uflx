"""MLIR code-generation backend for UFLx forms, built on the Python op-builder API.

See README.md for what this package does, its current restrictions, and
how to build the `mlir` Python bindings it depends on (not a pip
dependency -- see below).
"""

from __future__ import annotations

from uflx_mlir.emit import generate_mlir_module

__all__ = ["generate_mlir_module"]
