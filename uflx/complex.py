"""Complex number algorithms."""

from typing import Protocol, runtime_checkable

from uflx.algorithms import replace
from uflx.expressions import AbstractExpression, Conj
from uflx.graphs import GraphNode, as_graph


@runtime_checkable
class ComplexValued(Protocol):
    """A complex valued node."""

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""


def take_real_part(
    expression: GraphNode,
) -> GraphNode:
    """Take the real part of all complex values."""
    return replace(
        expression,
        {
            node: node.re
            for node in as_graph(expression)
            if isinstance(node, ComplexValued) and isinstance(node, GraphNode)
        },
    )


def take_imaginary_part(
    expression: GraphNode,
) -> GraphNode:
    """Take the imaginary part of all complex values."""
    return replace(
        expression,
        {
            node: node.im
            for node in as_graph(expression)
            if isinstance(node, ComplexValued) and isinstance(node, GraphNode)
        },
    )


def conj(value: AbstractExpression) -> AbstractExpression:
    """Get the complex conjugate."""
    if isinstance(value, ComplexValued):
        return value.re - value.im
    else:
        return Conj(value)
