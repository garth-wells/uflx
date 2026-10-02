# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from uflx.complex import conj
from uflx.domains import AbstractCoordinateElement, AbstractDomain
from uflx.expressions import AbstractExpression, BinaryOperator, UnaryOperator
from uflx.functions import AbstractFunction
from uflx.geometry import JacobianInverseTranspose
from uflx.graphs import GraphNode, as_graph
from uflx.maps import PushedForward
from uflx.tensors import Identity, Tensor, Vector, zero


class Inner(BinaryOperator):
    """Inner product operator.

    NOTE: document what happens here with conjugates.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class Transpose(UnaryOperator):
    """Transpose of a rank-2 (matrix-shaped) expression."""

    def __init__(self, argument: AbstractExpression):
        """Initialise."""
        if len(argument.value_shape) != 2:
            raise ValueError(
                "transpose() is only defined for rank-2 (matrix-shaped) expressions, got "
                f"shape {argument.value_shape}."
            )
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape[::-1]

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        i, j = indices
        return self.argument.component(j, i)


class Tr(UnaryOperator):
    """Trace of a square rank-2 (matrix-shaped) expression."""

    def __init__(self, argument: AbstractExpression):
        """Initialise."""
        shape = argument.value_shape
        if len(shape) != 2 or shape[0] != shape[1]:
            raise ValueError(
                "tr() is only defined for square rank-2 (matrix-shaped) expressions, got "
                f"shape {shape}."
            )
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return sum(
            self.argument.component(i, i).as_complex() for i in range(self.argument.value_shape[0])
        )

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return sum(
            self.argument.component(i, i).as_float() for i in range(self.argument.value_shape[0])
        )

    def as_int(self) -> int:
        """Convert to an integer."""
        return sum(
            self.argument.component(i, i).as_int() for i in range(self.argument.value_shape[0])
        )


class Grad(UnaryOperator):
    """Gradient operator."""

    def __init__(self, argument: GraphNode):
        """Initialise."""
        assert isinstance(argument, AbstractFunction) and not argument.is_reference
        self._physical_argument = argument
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        gdim = self._physical_argument.function_space.domain.geometric_dimension
        return (*self._physical_argument.value_shape, gdim)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError("Cannot get a 'component' of a Grad")

    def pull_back_to_reference(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the reference cell."""
        if self._physical_argument.is_cellwise_constant:
            # The gradient of a cellwise constant is exactly zero.
            return zero(self.value_shape)

        # assert isinstance(self.argument, EvaluatedPhysicalBasisFunction)
        argument = node_map.get(self.argument, self.argument)

        def extract_domain(node: GraphNode) -> AbstractDomain:
            """Extract the domain associated with a node."""
            domain: AbstractDomain | None = None
            for i in as_graph(node).descendants(node):
                if isinstance(i, AbstractFunction) and not i.is_reference:
                    if domain is None:
                        domain = i.function_space.domain
                    else:
                        assert domain == i.function_space.domain
            assert domain is not None
            return domain

        domain = extract_domain(self)
        assert isinstance(domain, AbstractCoordinateElement)
        if isinstance(argument, PushedForward):
            return JacobianInverseTranspose(domain) @ ReferenceGrad(argument.function)
        raise NotImplementedError()


class ReferenceGrad(UnaryOperator):
    """Gradient operator."""

    def __init__(self, argument: GraphNode):
        """Initialise."""
        assert isinstance(argument, AbstractFunction) and argument.is_reference
        self._reference_argument = argument
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (*self._reference_argument.value_shape, self._reference_argument.domain_size)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError(
            "Cannot get a 'component' of a ReferenceGrad. Try calling expand_geometry first"
        )

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        argument = self._reference_argument
        d = argument.domain_size
        if argument.is_cellwise_constant:
            return zero((*argument.value_shape, d))

        value_shape = argument.value_shape
        diffs = [argument.diff(i) for i in range(d)]
        if value_shape == ():
            # A scalar argument's directional derivatives are themselves
            # scalars -- nothing to index into.
            return Vector(diffs)

        # A non-scalar (eg vector-valued) argument's directional
        # derivatives are each still full value_shape-shaped expressions
        # (EvaluatedBasisFunction.diff keeps the same, un-indexed
        # component), so build the (*value_shape, d) result by indexing
        # into each direction's derivative rather than treating the
        # derivatives themselves as the leaves -- a plain
        # Vector(diffs)/Tensor(diffs) would silently misreport its own
        # shape as (d,), since Tensor's shape inference only looks at the
        # nesting of the Python list it is given, not each leaf's own
        # value_shape.
        def build(indices: tuple[int, ...]):
            if len(indices) == len(value_shape):
                return [diffs[i].component(*indices) for i in range(d)]
            return [build((*indices, k)) for k in range(value_shape[len(indices)])]

        return Tensor(build(()))


def grad(a: AbstractExpression) -> AbstractExpression:
    """The gradient of an expression."""
    if isinstance(a, AbstractFunction) and not a.is_reference and a.is_cellwise_constant:
        # The Grad of a cellwise constant physical function is zero.
        gdim = a.function_space.domain.geometric_dimension
        return zero((*a.value_shape, gdim))
    return Grad(a)


def inner(a: AbstractExpression, b: AbstractExpression) -> AbstractExpression:
    """Inner product."""
    if a.value_shape != b.value_shape:
        raise ValueError("Incompatible value shapes.")

    if a.value_shape == ():
        return a * conj(b)

    return Inner(a, b)


def transpose(a: AbstractExpression) -> AbstractExpression:
    """The transpose of a rank-2 (matrix-shaped) expression."""
    return Transpose(a)


def tr(a: AbstractExpression) -> AbstractExpression:
    """The trace of a square rank-2 (matrix-shaped) expression."""
    return Tr(a)


def sym(a: AbstractExpression) -> AbstractExpression:
    """The symmetric part of a square rank-2 (matrix-shaped) expression."""
    return (a + transpose(a)) / 2


def skew(a: AbstractExpression) -> AbstractExpression:
    """The skew-symmetric (antisymmetric) part of a square rank-2 expression."""
    return (a - transpose(a)) / 2


def dev(a: AbstractExpression) -> AbstractExpression:
    """The deviatoric (trace-free) part of a square rank-2 (matrix-shaped) expression."""
    d = a.value_shape[0]
    return a - (tr(a) / d) * Identity(d)
