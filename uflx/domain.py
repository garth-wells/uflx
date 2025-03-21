"""Domains."""

from abc import ABC, abstractmethod

from uflx.cell import AbstractCell


class AbstractDomain(ABC):
    """Abstract base class for a domain."""

    @property
    @abstractmethod
    def geometric_dimension(self) -> int:
        """The dimension of the space this domain is embedded in."""

    @property
    @abstractmethod
    def topological_dimension(self) -> int:
        """The dimension of the topology of this domain."""


class EmbeddedCell(AbstractDomain):
    """A cell embedded in R^d."""

    def __init__(self, cell: AbstractCell, gdim: int | None = None):
        """Initialise."""
        self._cell = cell
        if gdim is None:
            self._gdim = cell.topological_dimension
        else:
            self._gdim = gdim

    @property
    def geometric_dimension(self) -> int:
        return self._gdim

    @property
    def topological_dimension(self) -> int:
        return self._cell.topological_dimension
