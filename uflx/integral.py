"""Integral."""

from abc import ABC, abstractmethod

from uflx.expression import AbstractExpression


class AbstractMeasure(ABC):
    """Abstract base class for an integral measure."""

    def __rmul__(self, other):
        """Right multiply by an expression to form an integral."""
        if isinstance(other, AbstractExpression):
            return Integral(other, self)
        return NotImplemented


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


class Integral(AbstractIntegral):
    """An integral."""

    def __init__(self, integrand: AbstractExpression, measure: AbstractMeasure):
        """Initialise."""
        self._integrand = integrand
        self._measure = measure

    @property
    def integrand(self) -> AbstractExpression:
        """The integrand."""
        return self.integrand

    @property
    def measure(self) -> AbstractMeasure:
        """The integral measure."""
        return self._measure


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
