# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

"""Entities.

A entity is a single item in a mesh (usually a polytope). Entities of (topological)
dimension 0, 1, 2 and 3 are called vertices, edges, faces and volumes (respectively).

In a mesh, the entities of the highest dimension are called cells. The codimension
of an entity in this mesh is equal to the dimension of the entity subtracted from
the dimension of the cell. Entities of codimension 0, 1, 2 and 3 are called cells,
facets, ridges and peaks (respectively).
"""

from __future__ import annotations

from abc import ABC, abstractmethod


class AbstractEntity(ABC):
    """Abstract base class for entities."""

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check if this entity is equal to another entity."""

    @property
    @abstractmethod
    def topological_dimension(self) -> int:
        """Topological dimension of the entity."""

    @abstractmethod
    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            A list of sub-entities of the given dimension.
        """
