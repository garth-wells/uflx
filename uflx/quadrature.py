"""Quadrature rules."""

from typing import Self

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode
from uflx.utils import indented


class QuadratureRule:
    """A quadrature rule."""

    def __init__(self, points, weights):
        """Initialise."""
        self.points = points
        self.weights = weights

    @property
    def npoints(self):
        """The number of points in the quadrature rule."""
        return len(self.weights)


class QuadraturePoint:
    """A point in a quadrature rule."""

    def __init__(self, rule, index):
        """Initalise."""
        self.rule = rule
        self.index = index

    def __repr__(self):
        """Representation."""
        return f"QuadraturePoint({self.index})"


class QuadratureWeight(AbstractExpression):
    """A weight in a quadrature rule."""

    def __init__(self, rule, index):
        """Initalise."""
        self.rule = rule
        self.index = index

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def __repr__(self):
        """Representation."""
        return f"QuadratureWeight({self.index})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self.__class__(self.rule, self.index)


class QuadratureLoop:
    """A loop over the points in a quadrature rule."""

    def __init__(self, body, rule, variable):
        """Initalise."""
        self.body = body
        self.rule = rule
        self.variable = variable

    def __repr__(self):
        """Representation."""
        return f"QuadratureLoop({self.variable})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self.__class__(replacements.get(self.body, self.body), self.rule, self.variable)

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        return (
            f"for (int {self.variable}=0; {self.variable}!={self.rule.npoints}; "
            f"++{self.variable})\n"
            "{\n" + indented(self.body.generate_c(), 2) + "\n}"
        )
