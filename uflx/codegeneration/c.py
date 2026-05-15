"""Generation of C code."""

from typing import Protocol, runtime_checkable

import numpy as np


@runtime_checkable
class GenerateC(Protocol):
    """Protocol for Objects that can be converted to C code."""

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""


def c_table(table: np.ndarray) -> str:
    """Convert a numpy array to C."""
    if len(table.shape) == 1:
        return "{" + ", ".join(f"{i}" for i in table) + "}"
    return "{" + ", ".join(c_table(i) for i in table) + "}"


def tables_to_c(tables: dict[str, np.ndarray]) -> str:
    """Convert tables of values to a string of code."""
    return "\n".join(
        f"static const double {variable}["
        + "][".join(f"{i}" for i in table.shape)
        + "] = "
        + c_table(table)
        + ";"
        for variable, table in tables.items()
    )
