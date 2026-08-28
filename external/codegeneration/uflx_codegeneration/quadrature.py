"""Quadrature rules."""

from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode
from uflx.points import AbstractPoint, AbstractSetOfPoints

from uflx_codegeneration.c import GenerateC
from uflx_codegeneration.utils import indented


class QuadratureRule(AbstractSetOfPoints):
    """A quadrature rule."""

    def __init__(self, points: npt.NDArray[np.floating], weights: npt.NDArray[np.floating]):
        """Initialise."""
        assert points.shape[0] == len(weights)
        self.points = points
        self.weights = weights

    @property
    def npoints(self) -> int:
        """The number of points in the set."""
        return len(self.weights)

    @property
    def geometric_dimension(self) -> int:
        """The dimension of each point in the set."""
        return self.points.shape[1]


class QuadraturePoint(AbstractPoint):
    """A point in a quadrature rule."""

    def __init__(self, rule: QuadratureRule, index: int | str):
        """Initialise."""
        self.rule = rule
        self._index = index

    @property
    def points_set(self) -> QuadratureRule:
        """Get all the points in the set."""
        return self.rule

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return self.rule.geometric_dimension

    @property
    def index(self) -> int | str:
        """Get the index of the point in the set."""
        return self._index

    def __repr__(self):
        """Representation."""
        return f"QuadraturePoint({self._index})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

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

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


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

    def generate_c(self) -> str:
        """Generate code for this object."""
        if not isinstance(self.body, GenerateC):
            raise NotImplementedError(f"GenerateC is not implemented for {self.body.__class__}")
        return (
            f"for (int {self.variable}=0; {self.variable}!={self.rule.npoints}; "
            f"++{self.variable})\n"
            "{\n" + indented(self.body.generate_c(), 2) + "\n}"
        )


def quadrature_rule(
    points: Sequence[Sequence[float]] | npt.ArrayLike,
    weights: Sequence[float] | npt.ArrayLike,
) -> QuadratureRule:
    """Create a quadrature rule."""
    return QuadratureRule(np.array(points), np.array(weights))
