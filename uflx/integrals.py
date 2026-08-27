# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Measures and integrals."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import Graph, GraphNode, generate_graph


class AbstractMeasure(ABC):
    """Abstract base class for an integral measure."""

    def __rmul__(self, other: AbstractExpression) -> Integral:
        """Right multiply by an expression to form an integral."""
        if isinstance(other, AbstractExpression):
            return Integral(other, self)
        return NotImplemented

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""


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

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""

    def __eq__(self, other):
        """Check for equality."""
        if isinstance(other, AbstractIntegral):
            return self.integrand == other.integrand and self.measure == other.measure
        return NotImplemented

    def __hash__(self):
        """Hash."""
        return hash((hash(self.integrand), hash(self.measure)))

    def __repr__(self) -> str:
        """Representation."""
        return self.__class__.__name__


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

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._integrand, self._measure


class Measure(AbstractMeasure):
    """An integral measure."""

    def __init__(
        self, dim: int | None = None, codim: int | None = None, boundary_only: bool = False
    ):
        """Initialise."""
        self._dim = dim
        self._codim = codim
        self._boundary_only = boundary_only

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._dim, self._codim, self._boundary_only

    def __repr__(self) -> str:
        """Representation."""
        kwargs = {}
        if self._dim is not None:
            kwargs["dim"] = self._dim
        if self._codim is not None:
            kwargs["codim"] = self._codim
        if self._boundary_only:
            kwargs["boundary_only"] = self._boundary_only
        return (
            f"{self.__class__.__name__}("
            + ", ".join(f"{key}={value}" for key, value in kwargs.items())
            + ")"
        )


dx = Measure(codim=0)
