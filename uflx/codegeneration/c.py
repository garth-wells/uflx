"""Generation of C code."""

from typing import Protocol, Self, runtime_checkable


@runtime_checkable
class ConvertToCCode(Protocol):
    """Protocol for Objects that can be converted to C code."""

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
