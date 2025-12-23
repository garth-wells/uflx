"""Functions."""

from abc import abstractmethod

from uflx.expression import AbstractExpression
from uflx.function_space import AbstractFunctionSpace


class AbstractFunction(AbstractExpression):
    """Abstract base class for a function."""

    @property
    @abstractmethod
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""


class TestFunction(AbstractFunction):
    """A test function."""

    __test__ = False

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise."""
        self._space = space

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space


class TrialFunction(AbstractFunction):
    """A trial function."""

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise."""
        self._space = space

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space
