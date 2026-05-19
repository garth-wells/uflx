"""Generation of C code."""

from typing import Protocol, runtime_checkable

import numpy as np
from uflx.expressions import Abs, Add, Mult, Subtract
from uflx.quadrature import QuadratureLoop

from uflx_codegeneration.utils import indented


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


def mult_generate_c(self, bracketed: bool = False) -> str:
    """Generate code for this object."""
    assert isinstance(self.first, GenerateC)
    assert isinstance(self.second, GenerateC)
    return self.first.generate_c(True) + " * " + self.second.generate_c(True)


Mult.generate_c = mult_generate_c


def add_generate_c(self, bracketed: bool = False) -> str:
    """Generate code for this object."""
    assert isinstance(self.first, GenerateC)
    assert isinstance(self.second, GenerateC)
    code = self.first.generate_c() + " + " + self.second.generate_c()
    if bracketed:
        return f"({code})"
    else:
        return code


Add.generate_c = add_generate_c


def subtract_generate_c(self, bracketed: bool = False) -> str:
    """Generate code for this object."""
    assert isinstance(self.first, GenerateC)
    assert isinstance(self.second, GenerateC)
    code = self.first.generate_c() + " - " + self.second.generate_c()
    if bracketed:
        return f"({code})"
    else:
        return code


Subtract.generate_c = subtract_generate_c


def abs_generate_c(self, bracketed: bool = False) -> str:
    """Generate code for this object."""
    assert isinstance(self.argument, GenerateC)
    return f"fabs({self.argument.generate_c()})"


Abs.generate_c = abs_generate_c


def ql_generate_c(self, bracketed: bool = False) -> str:
    """Generate code for this object."""
    assert isinstance(self.body, GenerateC)
    return (
        f"for (int {self.variable}=0; {self.variable}!={self.rule.npoints}; "
        f"++{self.variable})\n"
        "{\n" + indented(self.body.generate_c(), 2) + "\n}"
    )


QuadratureLoop.generate_c = ql_generate_c
