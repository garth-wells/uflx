# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Measures and integrals."""

from abc import ABC, abstractmethod
from typing import Self

from uflx.expressions import AbstractExpression
from uflx.graphs import Graph, GraphNode, generate_graph


class AbstractMeasure(ABC):
    """Abstract base class for an integral measure."""

    def __rmul__(self, other):
        """Right multiply by an expression to form an integral."""
        if isinstance(other, AbstractExpression):
            return Integral(other, self)
        return NotImplemented

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self


class AbstractIntegral(ABC):
    """Abstract base class for an integral."""

    @property
    @abstractmethod
    def integrand(self) -> AbstractExpression:
        """The integrand."""

    @property
    @abstractmethod
    def measure(self) -> AbstractMeasure:
        """The integral measure."""

    @property
    def graph(self) -> Graph:
        """The graph that represents this object."""
        return generate_graph(self)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.integrand, self.measure}

    @abstractmethod
    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""


class Integral(AbstractIntegral):
    """An integral."""

    def __init__(self, integrand: AbstractExpression, measure: AbstractMeasure):
        """Initialise."""
        self._integrand = integrand
        self._measure = measure

    @property
    def integrand(self) -> AbstractExpression:
        """The integrand."""
        return self._integrand

    @property
    def measure(self) -> AbstractMeasure:
        """The integral measure."""
        return self._measure

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        if self._integrand not in replacements and self._measure not in replacements:
            return self
        return Integral(
            replacements.get(self._integrand, self._integrand),
            replacements.get(self._measure, self._measure),
        )


class Measure(AbstractMeasure):
    """An integral measure."""

    def __init__(self, *, dim: int | None = None, codim=int | None, boundary_only: bool = False):
        """Initialise."""
        self._dim = dim
        self._codim = codim
        self._boundary_only = boundary_only


dx = Measure(codim=0)
dS = Measure(codim=1)
ds = Measure(codim=1, boundary_only=True)
