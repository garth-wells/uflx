"""Quadrature rules."""

from abc import ABC, abstractmethod
from collections.abc import Hashable, Sequence
from typing import Self

import numpy as np

from uflx.codegeneration.c import GenerateC
from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode
from uflx.utils import indented


class PointInSet(ABC):
    """A single point in a set of points."""

    @property
    @abstractmethod
    def points(self) -> np.ndarray:
        """Get all the points in the set."""

    @property
    @abstractmethod
    def index(self) -> int | str:
        """Get the index of the point in the set."""

    @property
    @abstractmethod
    def points_id(self) -> Hashable:
        """Get an identifier for the set of points."""


class QuadratureRule:
    """A quadrature rule."""

    def __init__(self, points: np.ndarray, weights: np.ndarray):
        """Initialise."""
        self.points = points
        self.weights = weights

    @property
    def npoints(self):
        """The number of points in the quadrature rule."""
        return len(self.weights)


class QuadraturePoint(PointInSet):
    """A point in a quadrature rule."""

    def __init__(self, rule: QuadratureRule, index: int | str):
        """Initalise."""
        self.rule = rule
        self._index = index

    @property
    def index(self) -> int | str:
        """Get the index of the point in the set."""
        return self._index

    @property
    def points(self) -> np.ndarray:
        """Get all the points in the set."""
        return self.rule.points

    @property
    def points_id(self) -> Hashable:
        """Get an identifier for the set of points."""
        return self.rule

    def __repr__(self):
        """Representation."""
        return f"QuadraturePoint({self.index})"


class QuadratureWeight(AbstractExpression):
    """A weight in a quadrature rule."""

    def __init__(self, rule: QuadratureRule, index: int | str):
        """Initalise."""
        self.rule = rule
        self._index = index

    @property
    def index(self) -> int | str:
        """Get the index of the point in the set."""
        return self._index

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
        return self.__class__(self.rule, self._index)


class QuadratureLoop:
    """A loop over the points in a quadrature rule."""

    def __init__(self, body: GraphNode, rule: QuadratureRule, variable: str):
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
        assert isinstance(self.body, GenerateC)
        return (
            f"for (int {self.variable}=0; {self.variable}!={self.rule.npoints}; "
            f"++{self.variable})\n"
            "{\n" + indented(self.body.generate_c(), 2) + "\n}"
        )


def quadrature_rule(
    points: Sequence[Sequence[float]],
    weights: Sequence[float],
) -> QuadratureRule:
    """Create a quadrature rule."""
    return QuadratureRule(np.array(points), np.array(weights))
