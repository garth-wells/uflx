# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Finite element domains.

A domain is a subset of R^d over which something can be integrated.
"""

from abc import ABC, abstractmethod
from collections.abc import Sequence

from uflx.entities import AbstractEntity


class AbstractDomain(ABC):
    """Abstract base class for a domain."""

    @property
    @abstractmethod
    def geometric_dimension(self) -> int:
        """The dimension of the space this domain is embedded in."""

    @property
    @abstractmethod
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cells in the mesh."""


class Domain(AbstractDomain):
    """A set of cells embedded in R^d."""

    def __init__(self, cells: tuple[AbstractEntity, ...], gdim: int):
        """Initialise."""
        self._cells = cells
        self._gdim = gdim

    @property
    def geometric_dimension(self) -> int:
        """Dimension of the space this domain is embedded in."""
        return self._gdim

    @property
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cells in the mesh."""
        return self._cells


def domain(cells: Sequence[AbstractEntity] | AbstractEntity, gdim: int | None = None):
    """Create a domain.

    Args:
        cells: The cell or cells included in this domain
        gdim: The geometric dimension of the space in which this domain is embedded
    """
    if isinstance(cells, AbstractEntity):
        cells = (cells,)
    if gdim is None:
        gdim = max(cell.topological_dimension for cell in cells)

    return Domain(tuple(cells), gdim)
