"""Utility functions."""

from collections.abc import Sequence


def product(ls: Sequence[int]) -> int:
    """Return the product of numbers in a list."""
    result = 1
    for i in ls:
        result *= i
    return result
