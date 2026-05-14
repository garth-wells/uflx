"""Utility functions."""


def indented(code: str, spaces: int) -> str:
    """Add indentation to a block of code."""
    return "\n".join(" " * spaces + line for line in code.split("\n"))
