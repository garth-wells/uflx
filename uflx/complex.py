"""Complex number algorithms."""

from typing import Protocol, runtime_checkable
from uflx.graphs.algorithms import replace
from uflx.expressions import AbstractExpression


@runtime_checkable
class ComplexValued(Protocol):
    def re(self) -> AbstractExpression:
        """Get real part."""

    def im(self) -> AbstractExpression:
        """Get imaginary part."""


def take_real_part(
    graph: Graph,
) -> Graph:
    """Take the real part of all complex values."""
    return replace(graph, {node: node.re() for node in graph if isinstance(node, ComplexValued)})


def take_imaginary_part(
    graph: Graph,
) -> Graph:
    """Take the imaginary part of all complex values."""
    return replace(graph, {node: node.im() for node in graph if isinstance(node, ComplexValued)})
