# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Expression.

An expression is any algebraic expression that could be used as an integrand.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterable
from typing import Any

from uflx.graphs import GraphNode


class AbstractExpression(ABC):
    """Abstract base class for expressions."""

    @property
    @abstractmethod
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""

    @property
    @abstractmethod
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""

    def __mul__(self, other: Any) -> AbstractExpression:
        """Multiply."""
        if isinstance(other, AbstractExpression):
            if self.value_shape == other.value_shape:
                return Mult(self, other)
            if other.value_shape == ():
                return ScalarMult(other, self)
            if self.value_shape == ():
                return ScalarMult(self, other)
            raise ValueError(
                f"Cannot multiply expressions with shapes {self.value_shape} and "
                f"{other.value_shape}. To compute a matrix-vector or matrix-matrix"
                "product, use the '@' operator."
            )
        try:
            return to_scalar(other) * self
        except ValueError:
            return NotImplemented

    def __rmul__(self, other: Any) -> AbstractExpression:
        """Multiply."""
        try:
            return to_scalar(other) * self
        except ValueError:
            return NotImplemented

    def __truediv__(self, other: Any) -> AbstractExpression:
        """Division."""
        if isinstance(other, AbstractExpression):
            return Div(self, other)
        try:
            return self / to_scalar(other)
        except ValueError:
            return NotImplemented

    def __rtruediv__(self, other: Any) -> AbstractExpression:
        """Division."""
        try:
            return to_scalar(other) / self
        except ValueError:
            return NotImplemented

    def __add__(self, other: Any) -> AbstractExpression:
        """Add."""
        if isinstance(other, AbstractExpression):
            return Add(self, other)
        try:
            return self + to_scalar(other)
        except ValueError:
            return NotImplemented

    def __radd__(self, other: Any) -> AbstractExpression:
        """Add."""
        try:
            return to_scalar(other) + self
        except ValueError:
            return NotImplemented

    def __sub__(self, other: Any) -> AbstractExpression:
        """Subtract."""
        if isinstance(other, AbstractExpression):
            return Subtract(self, other)
        try:
            return self - to_scalar(other)
        except ValueError:
            return NotImplemented

    def __rsub__(self, other: Any) -> AbstractExpression:
        """Subtract."""
        try:
            return to_scalar(other) - self
        except ValueError:
            return NotImplemented

    def __neg__(self):
        """Negate."""
        return Neg(self)

    def __abs__(self):
        """Absolute value."""
        return Abs(self)

    def __repr__(self) -> str:
        """Representation."""
        return self.__class__.__name__

    def __pow__(self, power):
        """Raise to a power."""
        if isinstance(power, int):
            if power < 0:
                return RealScalar(1) / self**-power
            if power == 0:
                return RealScalar(1)
            return self * self ** (power - 1)
        return NotImplemented

    @abstractmethod
    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return Re(self)

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        return Im(self)

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        try:
            return complex(self.as_float())
        except ValueError:
            raise ValueError(f"Cannot convert {self.__class__.__name__} to complex")

    def as_float(self) -> float:
        """Convert to a floating point number."""
        try:
            return float(self.as_int())
        except ValueError:
            raise ValueError(f"Cannot convert {self.__class__.__name__} to float")

    def as_int(self) -> int:
        """Convert to an integer."""
        raise ValueError(f"Cannot convert {self.__class__.__name__} to int")


class AbstractScalar(AbstractExpression):
    """Abstract base class for scalars."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class AbstractInteger(AbstractScalar):
    """Abstract base class for integer values."""

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return self

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        return Integer(0)


class RealScalar(AbstractScalar):
    """A real scalar."""

    def __init__(self, value: float):
        """Initialise."""
        self.value = value

    def __repr__(self):
        """Representation."""
        return f"{self.value}"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.value,)

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return self

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        return RealScalar(0)

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return self.value


class ComplexScalar(AbstractScalar):
    """A complex scalar."""

    def __init__(self, re: AbstractScalar, im: AbstractScalar):
        """Initialise."""
        self._re = re
        self._im = im

    def __repr__(self):
        """Representation."""
        return f"{self._re!r} + ({self._im!r})j"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._re, self._im

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return self._re

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        return self._im

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return self._re.as_float() + 1j * self._im.as_float()


class Integer(AbstractInteger):
    """An integer."""

    def __init__(self, value: int):
        """Initialise."""
        self.value = value

    def __repr__(self):
        """Representation."""
        return f"{self.value}"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.value,)

    def as_int(self) -> int:
        """Convert to an integer."""
        return self.value


def to_scalar(value: Any) -> AbstractScalar:
    """Convert a value to a UFLx scalar or raise a ValueError if it cannot be converted."""
    if isinstance(value, float):
        return RealScalar(value)
    if isinstance(value, int):
        return Integer(value)
    if isinstance(value, complex):
        return ComplexScalar(RealScalar(value.real), RealScalar(value.imag))
    raise ValueError(f"Cannot convert value of type {type(value)} to UFLx scalar.")


class UnaryOperator(AbstractExpression):
    """A unary operator.

    Unary operators act on a single input.
    """

    def __init__(self, argument: AbstractExpression):
        """Initialise."""
        self.argument = argument

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.argument}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.argument,)

    def __repr__(self) -> str:
        """Representation."""
        return self.__class__.__name__


class Re(UnaryOperator):
    """Real part."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return Re(self.argument.component(*indices))

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return self

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        return RealScalar(0)

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return self.argument.as_complex().real


class Im(UnaryOperator):
    """Imaginary part."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return Im(self.argument.component(*indices))

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return RealScalar(0)

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        return self

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return self.argument.as_complex().imag


class BinaryOperator(AbstractExpression):
    """A binary operator.

    Binary operators act on two inputs.
    """

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        self.first = first
        self.second = second

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.first, self.second}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.first, self.second

    def __repr__(self) -> str:
        """Representation."""
        return self.__class__.__name__


class Mult(BinaryOperator):
    """Componentwise multiplication operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.first.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.first.component(*indices) * self.second.component(*indices)

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return self.first.as_complex() * self.second.as_complex()

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return self.first.as_float() * self.second.as_float()

    def as_int(self) -> int:
        """Convert to an integer."""
        return self.first.as_int() * self.second.as_int()


class ScalarMult(BinaryOperator):
    """Multiplication by a scalar."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.second.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.first * self.second.component(*indices)

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return self.first.as_complex() * self.second.as_complex()

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return self.first.as_float() * self.second.as_float()

    def as_int(self) -> int:
        """Convert to an integer."""
        return self.first.as_int() * self.second.as_int()


class Div(BinaryOperator):
    """Scalar multiplication operator."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        if second == 0:
            raise ZeroDivisionError()
        super().__init__(first, second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return self.first.as_complex() / self.second.as_complex()

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return self.first.as_float() / self.second.as_float()

    def as_int(self) -> int:
        """Convert to an integer."""
        return self.first.as_int() / self.second.as_int()


class Add(BinaryOperator):
    """Addition operator."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        assert first.value_shape == second.value_shape
        super().__init__(first, second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.first.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return self.first.component(*indices) + self.second.component(*indices)

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return self.first.re + self.second.re

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        return self.first.im + self.second.im

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return self.first.as_complex() + self.second.as_complex()

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return self.first.as_float() + self.second.as_float()

    def as_int(self) -> int:
        """Convert to an integer."""
        return self.first.as_int() + self.second.as_int()


class Subtract(BinaryOperator):
    """Subtraction operator."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        assert first.value_shape == second.value_shape
        super().__init__(first, second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.first.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return self.first.component(*indices) - self.second.component(*indices)

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return self.first.as_complex() - self.second.as_complex()

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return self.first.as_float() - self.second.as_float()

    def as_int(self) -> int:
        """Convert to an integer."""
        return self.first.as_int() - self.second.as_int()


class Abs(UnaryOperator):
    """Absolute value operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return Abs(self.argument.component(*indices))

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return abs(self.argument.as_float())

    def as_int(self) -> int:
        """Convert to an integer."""
        return abs(self.argument.as_int())


class Neg(UnaryOperator):
    """Negation operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return Neg(self.argument.component(*indices))

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return -self.argument.as_complex()

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return -self.argument.as_float()

    def as_int(self) -> int:
        """Convert to an integer."""
        return -self.argument.as_int()


class MatVec(BinaryOperator):
    """Matrix-vector multiplication operator."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        assert (
            len(first.value_shape) == 2
            and len(second.value_shape) == 1
            and first.value_shape[0] == second.value_shape[0]
        )
        super().__init__(first, second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.second.value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        (index,) = indices
        return expression_sum(
            self.first.component(index, i) * self.second.component(i)
            for i in range(self.first.value_shape[1])
        )


def expression_sum(
    expressions: Iterable[AbstractExpression], default: AbstractExpression | None = None
):
    """Take the sum of a sequence of expressions."""
    result = None
    for e in expressions:
        if result is None:
            result = e
        else:
            result += e
    if result is None:
        if default is None:
            raise ValueError("Cannot sum an empty sequence without a default return value")
        return default
    return result
