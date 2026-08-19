"""Quadrature rules."""

from collections.abc import Hashable, Sequence
from typing import Any

import numpy as np

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode
from uflx.points import PointInSet


class QuadratureRule:
    """A quadrature rule."""

    def __init__(self, points: npt.NDArray[np.floating], weights: npt.NDArray[np.floating]):
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
        """Initialise."""
        self.rule = rule
        self._index = index

    @property
    def index(self) -> int | str:
        """Get the index of the point in the set."""
        return self._index

    @property
    def points(self) -> npt.NDArray[np.floating]:
        """Get all the points in the set."""
        return self.rule.points

    @property
    def points_id(self) -> Hashable:
        """Get an identifier for the set of points."""
        return self.rule

    def __repr__(self):
        """Representation."""
        return f"QuadraturePoint({self._index})"

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.rule, self._index


class QuadratureWeight(AbstractExpression):
    """A weight in a quadrature rule."""

    def __init__(self, rule: QuadratureRule, index: int | str):
        """Initialise."""
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

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.rule, self._index


class QuadratureLoop:
    """A loop over the points in a quadrature rule."""

    def __init__(self, body: GraphNode, rule: QuadratureRule, variable: str):
        """Initialise."""
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

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.body, self.rule, self.variable


def quadrature_rule(
    points: Sequence[Sequence[float]],
    weights: Sequence[float],
) -> QuadratureRule:
    """Create a quadrature rule."""
    return QuadratureRule(np.array(points), np.array(weights))
