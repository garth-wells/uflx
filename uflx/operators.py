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
from uflx.tensors import Vector, zero


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

    def expand_geometry(self) -> GraphNode:
        """Expand geometry."""
        argument = self._reference_argument
        if argument.is_cellwise_constant:
            return zero((*argument.value_shape, argument.domain_size))
        return Vector([argument.diff(i) for i in range(argument.domain_size)])


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
