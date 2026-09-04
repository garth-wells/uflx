# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from uflx.expressions import AbstractExpression, BinaryOperator, UnaryOperator
from uflx.functions import AbstractPhysicalFunction, AbstractReferenceFunction
from uflx.geometry import JacobianInverseTranspose
from uflx.graphs import GraphNode, generate_graph
from uflx.tensors import Vector


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
        assert isinstance(argument, AbstractPhysicalFunction)
        self._physical_argument = argument
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._physical_argument.function_space.domain.geometric_dimension,)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError("Cannot get a 'component' of a Grad")

    def pull_back_to_reference(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the reference cell."""
        # assert isinstance(self.argument, EvaluatedPhysicalBasisFunction)
        argument = node_map.get(self.argument, self.argument)

        def extract_domain(node: GraphNode) -> AbstractDomain:
            """Extract the domain associated with a node."""
            domain: AbstractDomain | None = None
            for i in generate_graph(node).descendants(node):
                if isinstance(i, AbstractPhysicalFunction):
                    if domain is None:
                        domain = i.function_space.domain
                    else:
                        assert domain == i.function_space.domain
            assert domain is not None
            return domain

        domain = extract_domain(self)
        # assert isinstance(domain, AbstractCoordinateElement)
        # if isinstance(argument, PushedForward):
        return JacobianInverseTranspose(
            domain,
            None,
        ) * ReferenceGrad(argument.function)


class ReferenceGrad(UnaryOperator):
    """Gradient operator."""

    def __init__(self, argument: GraphNode):
        """Initialise."""
        assert isinstance(argument, AbstractReferenceFunction)
        self._reference_argument = argument
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._reference_argument.domain_size,)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError(
            "Cannot get a 'component' of a ReferenceGrad. Try calling expand_geometry first"
        )

    def expand_geometry(self) -> GraphNode:
        """Expand geometry."""
        argument = self._reference_argument
        return Vector([argument.diff(i) for i in range(argument.domain_size)])


class Conj(UnaryOperator):
    """Complex conjugate operator."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape

    def re(self) -> AbstractExpression:
        """Get real part."""
        return self.argument

    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        raise NotImplementedError()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        if self.value_shape == ():
            raise NotImplementedError("Cannot get a 'component' of a Grad")
        return Conj(self.argument.component(*indices))


def grad(a: AbstractExpression) -> Grad:
    """The gradient of an expression."""
    return Grad(a)


def inner(a: AbstractExpression, b: AbstractExpression) -> AbstractExpression:
    """Inner product."""
    if a.value_shape != b.value_shape:
        raise ValueError("Incompatible value shapes.")

    if a.value_shape == ():
        return a * Conj(b)

    return Inner(a, b)
