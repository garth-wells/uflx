"""Generation of C code."""

from typing import Protocol, runtime_checkable

import numpy as np
from uflx.expressions import Abs, Add, Mult, Subtract
from uflx.geometry import CoordinateDofComponent
from uflx.points import PointComponent

from uflx_codegeneration import symbols


@runtime_checkable
class GenerateC(Protocol):
    """Protocol for Objects that can be converted to C code."""

    def generate_c(self) -> str:
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


def mult_generate_c(self) -> str:
    """Generate code for this object."""
    if not isinstance(self.first, GenerateC):
        raise NotImplementedError(f"GenerateC is not implemented for {self.first.__class__}")
    if not isinstance(self.second, GenerateC):
        raise NotImplementedError(f"GenerateC is not implemented for {self.second.__class__}")
    return f"({self.first.generate_c()} * {self.second.generate_c()})"


setattr(Mult, "generate_c", mult_generate_c)


def add_generate_c(self) -> str:
    """Generate code for this object."""
    if not isinstance(self.first, GenerateC):
        raise NotImplementedError(f"GenerateC is not implemented for {self.first.__class__}")
    if not isinstance(self.second, GenerateC):
        raise NotImplementedError(f"GenerateC is not implemented for {self.second.__class__}")
    return f"({self.first.generate_c()} + {self.second.generate_c()})"


setattr(Add, "generate_c", add_generate_c)


def subtract_generate_c(self) -> str:
    """Generate code for this object."""
    if not isinstance(self.first, GenerateC):
        raise NotImplementedError(f"GenerateC is not implemented for {self.first.__class__}")
    if not isinstance(self.second, GenerateC):
        raise NotImplementedError(f"GenerateC is not implemented for {self.second.__class__}")
    return f"({self.first.generate_c()} - {self.second.generate_c()})"


setattr(Subtract, "generate_c", subtract_generate_c)


def abs_generate_c(self) -> str:
    """Generate code for this object."""
    if not isinstance(self.argument, GenerateC):
        raise NotImplementedError(f"GenerateC is not implemented for {self.argument.__class__}")
    return f"fabs({self.argument.generate_c()})"


setattr(Abs, "generate_c", abs_generate_c)


def pc_generate_c(self) -> str:
    """Generate code for this object."""
    c = self.point.get_component(self.component)
    if isinstance(c, PointComponent):
        raise NotImplementedError(f"GenerateC is not implemented for {self.__class__}")
    if not isinstance(c, GenerateC):
        raise NotImplementedError(f"GenerateC is not implemented for {c.__class__}")
    return c.generate_c()


setattr(PointComponent, "generate_c", pc_generate_c)


def cdc_generate_c(self) -> str:
    """Generate code for this object."""
    return f"{symbols.coordinate_dofs}[{self._tdim * self._point + self._component}]"


setattr(CoordinateDofComponent, "generate_c", cdc_generate_c)
