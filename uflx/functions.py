# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Functions.

A function is an item contained in a function space.
"""

from __future__ import annotations

from abc import abstractmethod
from itertools import count
from typing import Any, Self

from uflx.expressions import AbstractExpression
from uflx.function_spaces import AbstractFunctionSpace, AbstractReferenceMappedFunctionSpace
from uflx.graphs import GraphNode
from uflx.maps import PushedForward


class AbstractFunction(AbstractExpression):
    """Abstract base class for a function."""

    @property
    @abstractmethod
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""

    @abstractmethod
    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""


class AbstractPhysicalFunction(AbstractFunction):
    """Abstract base class for a function on a physical cell."""

    @property
    @abstractmethod
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.function_space.value_shape

    @property
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""
        return self.function_space.domain.cells[0].topological_dimension


class AbstractReferenceFunction(AbstractFunction):
    """Abstract base class for a function on a reference cell."""

    @property
    @abstractmethod
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""


class AbstractIntegralScopedFunction(AbstractPhysicalFunction):
    """Base class for a physical function that gets relabelled per integral.

    Shared by Argument and Coefficient. When an expression is wrapped in an
    Integral (see Integral.__init__), every not-yet-labelled instance of this
    class that the expression contains is reconstructed with that integral's
    label, so that -- if the same TestFunction/TrialFunction/Coefficient
    object is reused across more than one integral in a form -- each integral
    gets its own independently-labelled copy, and code generation (which
    looks up "arguments/coefficients whose integral_label matches this
    integral") never conflates the two.
    """

    def __init__(self, space: AbstractFunctionSpace, integral_label: str | None = None):
        """Initialise."""
        self._space = space
        self._integral_label = integral_label

    @property
    def integral_label(self) -> str | None:
        """Get the label of the integral that this function is included in."""
        return self._integral_label

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @abstractmethod
    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct this function with the given integral label."""


class AbstractReferenceIntegralScopedFunction(AbstractReferenceFunction):
    """Base class for a reference function that gets relabelled per integral."""

    def __init__(self, space: AbstractFunctionSpace, integral_label: str | None = None):
        """Initialise."""
        self._space = space
        self._integral_label = integral_label

    @property
    def integral_label(self) -> str | None:
        """Get the label of the integral that this function is included in."""
        return self._integral_label

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @abstractmethod
    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct this function with the given integral label."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        assert isinstance(self.function_space, AbstractReferenceMappedFunctionSpace)
        return self.function_space.elements[0].reference_value_shape


class Argument(AbstractIntegralScopedFunction):
    """A function that is a dimension of the tensor to be assembled."""

    def __init__(
        self, space: AbstractFunctionSpace, component: int, integral_label: str | None = None
    ):
        """Initialise.

        Args:
            space: The function space that this function lives in
            component: The component of the finite element tensor
                       to be assembled that this function represents
            integral_label: The label of the integral that this
                            argument is included in
        """
        super().__init__(space, integral_label)
        self._component = component

    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct the argument with the given integral label."""
        return self.__class__(self._space, self._component, integral_label)

    @property
    def component_index(self) -> int:
        """The component of the finite element tensor that this function represents."""
        return self._component

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self._component, self.integral_label

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError()

    def get_replacement(self, replacements: dict[GraphNode, GraphNode]) -> GraphNode | None:
        """Get the node to replace this node with, or None if no replacement can be made."""
        for old, new in replacements.items():
            if (
                isinstance(old, Argument)
                and old.function_space == self.function_space
                and old.component_index == self.component_index
            ):
                if (
                    isinstance(new, Argument)
                    and self.integral_label is not None
                    and new.integral_label is None
                ):
                    return new.reconstruct_with_integral_label(self.integral_label)
                return new


class ReferenceArgument(AbstractReferenceIntegralScopedFunction):
    """A function that is a dimension of the tensor to be assembled on the reference cell."""

    def __init__(
        self, space: AbstractFunctionSpace, component: int, integral_label: str | None = None
    ):
        """Initialise.

        Args:
            space: The function space that this function lives in
            component: The component of the finite element tensor
                       to be assembled that this function represents
            integral_label: The label of the integral that this
                            argument is included in
        """
        self._space = space
        self._component = component
        self._integral_label = integral_label

    @property
    def integral_label(self) -> str | None:
        """Get the label of the integral that this argument is included in."""
        return self._integral_label

    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct the argument with the given integral label."""
        return self.__class__(self._space, self._component, integral_label)

    @property
    def component_index(self) -> int:
        """The component of the finite element tensor that this function represents."""
        return self._component

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self._component

    @property
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""
        return self._space.domain.cells[0].topological_dimension

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError()


class TestFunction(Argument):
    """A test function."""

    __test__ = False

    def __init__(self, space: AbstractFunctionSpace, integral_label: str | None = None):
        """Initialise."""
        super().__init__(space, 0, integral_label)

    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct the argument with the given integral label."""
        return self.__class__(self._space, integral_label)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self.integral_label

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()

    def pull_back_to_reference(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the reference cell."""
        assert isinstance(self._space, AbstractReferenceMappedFunctionSpace)
        return PushedForward(
            self._space.elements[0].reference_map, ReferenceTestFunction(self._space)
        )


class TrialFunction(Argument):
    """A trial function."""

    def __init__(self, space: AbstractFunctionSpace, integral_label: str | None = None):
        """Initialise."""
        super().__init__(space, 1, integral_label)

    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct the argument with the given integral label."""
        return self.__class__(self._space, integral_label)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self.integral_label

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()

    def pull_back_to_reference(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the reference cell."""
        assert isinstance(self._space, AbstractReferenceMappedFunctionSpace)
        return PushedForward(
            self._space.elements[0].reference_map, ReferenceTrialFunction(self._space)
        )


class Coefficient(AbstractIntegralScopedFunction):
    """A known function with given degree-of-freedom values.

    Unlike Argument (a bound variable of the bilinear/linear form being
    assembled, which becomes an axis of the assembled tensor), a Coefficient
    represents an already-known function -- eg a previous solution, a
    material property, or any other field supplied at assembly time -- whose
    value is fully determined by a fixed array of degree-of-freedom values,
    not by the tensor being assembled.
    """

    _n = count(0)

    def __init__(self, space: AbstractFunctionSpace):
        """Initialise.

        Args:
            space: The function space that this function lives in
        """
        super().__init__(space)
        self._count = next(self._n)

    @classmethod
    def _restore(
        cls,
        space: AbstractFunctionSpace,
        count: int,
        integral_label: str | None,
    ) -> Self:
        """Restore a coefficient without allocating a new identity."""
        coefficient = cls.__new__(cls)
        AbstractIntegralScopedFunction.__init__(coefficient, space, integral_label)
        coefficient._count = count
        return coefficient

    @property
    def count(self) -> int:
        """A value that, together with function_space, uniquely identifies this coefficient."""
        return self._count

    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct the coefficient with the given integral label."""
        return self.__class__._restore(self._space, self._count, integral_label)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._space,)

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError()

    def get_replacement(self, replacements: dict[GraphNode, GraphNode]) -> GraphNode | None:
        """Get the node to replace this node with, or None if no replacement can be made."""
        for old, new in replacements.items():
            if (
                isinstance(old, Coefficient)
                and old.function_space == self.function_space
                and old.count == self.count
            ):
                if (
                    isinstance(new, Coefficient)
                    and self.integral_label is not None
                    and new.integral_label is None
                ):
                    return new.reconstruct_with_integral_label(self.integral_label)
                return new


class ReferenceTestFunction(ReferenceArgument):
    """A test function on the reference cell."""

    def __init__(self, space: AbstractFunctionSpace, integral_label: str | None = None):
        """Initialise."""
        super().__init__(space, 0, integral_label)

    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct the argument with the given integral label."""
        return self.__class__(self._space, integral_label)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._space, self._integral_label)

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()


class ReferenceTrialFunction(ReferenceArgument):
    """A trial function on the reference cell."""

    def __init__(self, space: AbstractFunctionSpace, integral_label: str | None = None):
        """Initialise."""
        super().__init__(space, 1, integral_label)

    def reconstruct_with_integral_label(self, integral_label: str) -> Self:
        """Reconstruct the argument with the given integral label."""
        return self.__class__(self._space, integral_label)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._space,)

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()
