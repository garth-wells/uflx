"""Cell."""

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

    @property
    def vertices(self) -> list[AbstractCell]:
        """Get the vertices of the cell."""
        return self.sub_entities(0)

    @property
    def edges(self) -> list[AbstractCell]:
        """Get the edges of the cell."""
        return self.sub_entities(1)

    @property
    def faces(self) -> list[AbstractCell]:
        """Get the faces of the cell."""
        return self.sub_entities(2)

    @property
    def facets(self) -> list[AbstractCell]:
        """Get the facets of the cell.

        Facets are the sub-entities of dimension tdim - 1.
        """
        return self.sub_entities(self.topological_dimension - 1)

    @property
    def ridges(self) -> list[AbstractCell]:
        """Get the ridges of the cell.

        Ridges are the sub-entities of dimension tdim - 2.
        """
        return self.sub_entities(self.topological_dimension - 2)

    @property
    def peaks(self) -> list[AbstractCell]:
        """Get the peaks of the cell.

        Peaks are the sub-entities of dimension tdim - 3.
        """
        return self.sub_entities(self.topological_dimension - 3)
