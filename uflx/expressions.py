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
from collections.abc import Iterable, Sequence
from math import gcd, prod
from typing import Any

from uflx.graphs.graphs import GraphNode


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
                return Product([self, other])
            if other.value_shape == ():
                return ScalarMult(other, self)
            if self.value_shape == ():
                return ScalarMult(self, other)
            raise ValueError(
                f"Cannot multiply expressions with shapes {self.value_shape} and "
                f"{other.value_shape}. To compute a matrix-vector or matrix-matrix "
                "product, use the '@' operator."
            )
        try:
            return to_scalar(other) * self
        except ValueError:
            return NotImplemented

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, self.__class__) and self.init_args == other.init_args

    def __hash__(self) -> int:
        """Hash."""
        return hash((f"uflx.{self.__class__.__name__}", *self.init_args))

    def __matmul__(self, other: Any) -> AbstractExpression:
        """Matrix multiply."""
        if isinstance(other, AbstractExpression):
            if self.value_shape[-1] != other.value_shape[0]:
                raise ValueError("Incompatible dimensions in matmul.")
            return MatMult(self, other)
        return NotImplemented

    def __rmul__(self, other: Any) -> AbstractExpression:
        """Multiply."""
        try:
            return to_scalar(other) * self
        except ValueError:
            return NotImplemented

    def __recip__(self) -> AbstractExpression:
        """Reciprocal."""
        return Reciprocal(self)

    def __truediv__(self, other: Any) -> AbstractExpression:
        """Division."""
        if isinstance(other, AbstractExpression):
            return self * other.__recip__()
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
            return Sum([self, other])
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
            return self + -other
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

    def __neg__(self) -> AbstractExpression:
        """Negate."""
        return Neg(self)

    def __abs__(self) -> AbstractExpression:
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

    def __recip__(self) -> AbstractExpression:
        """Reciprocal."""
        return RealScalar(1 / self.value)

    def __neg__(self) -> AbstractExpression:
        """Negation."""
        return RealScalar(-self.value)

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

    def __init__(self, real_part: AbstractScalar, imag_part: AbstractScalar):
        """Initialise."""
        self._re = real_part
        self._im = imag_part

    def __repr__(self):
        """Representation."""
        return f"{self._re!r} + ({self._im!r})j"

    def __recip__(self) -> AbstractExpression:
        """Reciprocal."""
        n = self._re**2 + self._im**2
        re = self._re / n
        im = self._im / n
        if isinstance(re, AbstractScalar) and isinstance(im, AbstractScalar):
            return ComplexScalar(re, im)
        else:
            return Div(Integer(1), self)

    def __neg__(self) -> AbstractExpression:
        """Negation."""
        re = -self.re
        im = -self.im
        if isinstance(re, AbstractScalar) and isinstance(im, AbstractScalar):
            return ComplexScalar(re, im)
        else:
            return Neg(self)

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

    def simplified_product(self, other: AbstractExpression) -> AbstractExpression | None:
        """Return a single expression representing the simplified product.

        This function should return None if no simplification can be made.
        """
        from uflx.functions import Argument

        if self.value == 1:
            return other
        if self.value == 0 and not isinstance(other, Argument):
            return self
        if isinstance(other, Integer):
            return Integer(self.value * other.value)

    def simplified_sum(self, other: AbstractExpression) -> AbstractExpression | None:
        """Return a single expression representing the simplified sum.

        This function should return None if no simplification can be made.
        """
        if self.value == 0:
            return other
        if isinstance(other, Integer):
            return Integer(self.value + other.value)

    def __eq__(self, other):
        """Check for equality."""
        if isinstance(other, Integer):
            return self.value == other.value
        return self.value == other

    def __hash__(self):
        """Hash."""
        return hash(("uflx.Integer", self.value))

    def __repr__(self):
        """Representation."""
        return f"{self.value}"

    def __recip__(self) -> AbstractExpression:
        """Reciprocal."""
        return Rational(1, self.value)

    def __neg__(self) -> AbstractExpression:
        """Negation."""
        return Integer(-self.value)

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


class Rational(AbstractScalar):
    """A rational number."""

    def __init__(self, numerator: int, denominator: int):
        """Initialise."""
        factor = gcd(numerator, denominator)
        self.numerator = numerator // factor
        self.denominator = denominator // factor

    def simplified_product(self, other: AbstractExpression) -> AbstractExpression | None:
        """Return a single expression representing the simplified product.

        This function should return None if no simplification can be made.
        """
        if isinstance(other, Integer):
            numerator = self.numerator * other.value
            denominator = self.denominator
        elif isinstance(other, Rational):
            numerator = self.numerator * other.numerator
            denominator = self.denominator * other.denominator
        else:
            return None

        factor = gcd(numerator, denominator)
        numerator //= factor
        denominator //= factor

        if denominator == 1 or numerator == 0:
            return Integer(numerator)
        else:
            return Rational(numerator, denominator)

    def simplified_sum(self, other: AbstractExpression) -> AbstractExpression | None:
        """Return a single expression representing the simplified product.

        This function should return None if no simplification can be made.
        """
        if isinstance(other, Integer):
            numerator = self.numerator + other.value * self.denominator
            denominator = self.denominator
        elif isinstance(other, Rational):
            numerator = self.numerator * other.denominator + self.denominator * other.numerator
            denominator = self.denominator * other.denominator
        else:
            return None

        factor = gcd(numerator, denominator)
        numerator //= factor
        denominator //= factor

        if denominator == 1 or numerator == 0:
            return Integer(numerator)
        else:
            return Rational(numerator, denominator)

    def __eq__(self, other):
        """Check for equality."""
        if isinstance(other, Rational):
            return self.numerator == other.numerator and self.denominator == other.denominator
        return False

    def __hash__(self):
        """Hash."""
        return hash(("uflx.Rational", self.numerator, self.denominator))

    def __repr__(self):
        """Representation."""
        return f"{self.numerator}/{self.denominator}"

    def __recip__(self) -> AbstractExpression:
        """Reciprocal."""
        return Rational(self.denominator, self.numerator)

    def __neg__(self) -> AbstractExpression:
        """Negation."""
        return Rational(-self.numerator, self.denominator)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.numerator, self.denominator)

    def as_float(self) -> float:
        """Convert to a float."""
        return self.numerator / self.denominator


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


class Conj(UnaryOperator):
    """Complex conjugate operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        return self.argument

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        raise NotImplementedError()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise NotImplementedError("Cannot get a 'component' of a Grad")
        return Conj(self.argument.component(*indices))


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


class Product(AbstractExpression):
    """Componentwise product."""

    def __init__(self, items: Sequence[AbstractExpression]):
        """Initialise."""
        if len(items) == 0:
            raise ValueError("Cannot create an empty product")
        self._value_shape = items[0].value_shape
        for i in items[1:]:
            assert i.value_shape == self._value_shape
        self._items = tuple(items)

    def __mul__(self, other: Any) -> AbstractExpression:
        """Multiply."""
        if isinstance(other, AbstractExpression):
            return Product(self._items + (other._items if isinstance(other, Product) else (other,)))
        try:
            return self * to_scalar(other)
        except ValueError:
            return NotImplemented

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set(self._items)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._items,)

    def __repr__(self) -> str:
        """Representation."""
        return "Product([" + ", ".join(f"{i!r}" for i in self._items) + "])"

    def __str__(self) -> str:
        """Representation."""
        return "Product"

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return Product([i.component(*indices) for i in self._items])

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return prod(i.as_complex() for i in self._items)

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return prod(i.as_float() for i in self._items)

    def as_int(self) -> int:
        """Convert to an integer."""
        return prod(i.as_int() for i in self._items)


class Sum(AbstractExpression):
    """Componentwise sum."""

    def __init__(self, items: Sequence[AbstractExpression]):
        """Initialise."""
        if len(items) == 0:
            raise ValueError("Cannot create an empty sum")
        self._value_shape = items[0].value_shape
        for i in items[1:]:
            assert i.value_shape == self._value_shape
        self._items = tuple(items)

    def __add__(self, other: Any) -> AbstractExpression:
        """Add."""
        if isinstance(other, AbstractExpression):
            return Sum(self._items + (other._items if isinstance(other, Sum) else (other,)))
        try:
            return self + to_scalar(other)
        except ValueError:
            return NotImplemented

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set(self._items)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._items,)

    def __repr__(self) -> str:
        """Representation."""
        return "Sum([" + ", ".join(f"{i!r}" for i in self._items) + "])"

    def __str__(self) -> str:
        """Representation."""
        return "Sum"

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._value_shape

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return Sum([i.component(*indices) for i in self._items])

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return sum(i.as_complex() for i in self._items)

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return sum(i.as_float() for i in self._items)

    def as_int(self) -> int:
        """Convert to an integer."""
        return sum(i.as_int() for i in self._items)


class ScalarMult(BinaryOperator):
    """Multiplication by a scalar."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        assert first.value_shape == ()
        assert second.value_shape != ()
        super().__init__(first, second)

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


class MatMult(BinaryOperator):
    """Multiplication by a matrix."""

    def __init__(self, first: AbstractExpression, second: AbstractExpression):
        """Initialise."""
        assert first.value_shape[-1] == second.value_shape[0]
        super().__init__(first, second)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.first.value_shape[:-1] + self.second.value_shape[1:]

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        assert len(indices) == len(self.value_shape)
        n = len(self.first.value_shape) - 1
        return expression_sum(
            self.first.component(*indices[:n], i) * self.second.component(i, *indices[n:])
            for i in range(self.first.value_shape[-1])
        )


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
        a = self.first.as_int()
        b = self.second.as_int()
        if a % b != 0:
            raise ValueError(f"Cannot convert {self.__class__.__name__} to int")
        return a // b


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
        try:
            return abs(self.argument.as_float())
        except ValueError:
            return abs(self.argument.as_complex())

    def as_int(self) -> int:
        """Convert to an integer."""
        return abs(self.argument.as_int())


class Reciprocal(UnaryOperator):
    """Reciprocal operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def simplified_product(self, other: AbstractExpression) -> AbstractExpression | None:
        """Return a single expression representing the simplified product.

        This function should return None if no simplification can be made.
        """
        if self.value_shape == () and other == self.argument:
            return Integer(1)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return Reciprocal(self.argument.component(*indices))

    def as_complex(self) -> complex:
        """Convert to a floating point number."""
        return 1.0 / self.argument.as_complex()

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return 1.0 / self.argument.as_float()

    def as_int(self) -> int:
        """Convert to an integer."""
        i = self.argument.as_int()
        if i in [-1, 1]:
            return i
        raise ValueError(f"Cannot convert {self.__class__.__name__} to int")


class Neg(UnaryOperator):
    """Negation operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def simplified_sum(self, other: AbstractExpression) -> AbstractExpression | None:
        """Return a single expression representing the simplified sum.

        This function should return None if no simplification can be made.
        """
        if self.value_shape == () and other == self.argument:
            return Integer(0)

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
