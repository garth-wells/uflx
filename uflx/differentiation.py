"""Differentitation of expressions."""

from typing import Protocol, runtime_checkable

from uflx.operators import UnaryOperator
from uflx.expressions import AbstractExpression


@runtime_checkable
class CanGrad(Protocol):
    """An expression whose grad can be done explicitly."""

    def grad(self) -> AbstractExpression:
        """The gradient."""



class Grad(UnaryOperator):
    """Gradient operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.argument.function_space.domain.geometric_dimension,)  # type: ignore


def grad(a: AbstractExpression) -> AbstractExpression:
    """The gradient of an expression."""
    if isinstance(a, CanGrad):
        return a.grad()

    return Grad(a)
