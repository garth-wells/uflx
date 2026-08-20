"""Utility functions."""

from uflx.entities import AbstractEntity


def indented(code: str, spaces: int) -> str:
    """Add indentation to a block of code."""
    return "\n".join(" " * spaces + line for line in code.split("\n"))


def index(*indices: int):
    """Get the index from a set of indices using triangular (2D) or tetrahedral (3D) ordering."""
    match len(indices):
        case 0:
            return 0
        case 1:
            (i,) = indices
            return i
        case 2:
            i, j = indices
            return (i + j + 1) * (i + j) // 2 + j
        case 3:
            i, j, k = indices
            return (
                (i + j + k) * (i + j + k + 1) * (i + j + k + 2) // 6
                + (j + k) * (j + k + 1) // 2
                + k
            )
        case _:
            raise ValueError("Invalid number of indices.")


def number_of_derivatives(nderivs: int, cell: AbstractEntity):
    """Get the total number of partial derivatives of degree <= nderivs."""
    return index(nderivs + 1, *[0 for i in range(1, cell.topological_dimension)])
