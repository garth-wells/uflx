# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

"""Finite element cells Cell."""

from __future__ import annotations

from abc import ABC, abstractmethod


class AbstractCell(ABC):
    """Abstract base class for cells."""

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check if this cell is equal to another cell."""

    @property
    @abstractmethod
    def topological_dimension(self) -> int:
        """The topological dimension of the cell."""

    @abstractmethod
    def sub_entities(self, dim: int) -> list[AbstractCell]:
        """Get a list of sub-entities of a given dimension."""
