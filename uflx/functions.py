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
from typing import Any, Self, cast

from uflx.domains import AbstractDomain, AbstractFiniteElementDomain
from uflx.expressions import AbstractExpression, Im, Re
from uflx.function_spaces import AbstractFunctionSpace, AbstractReferenceMappedFunctionSpace
from uflx.graphs import GraphNode
from uflx.maps import PushedForward
from uflx.tensors import zero


class AbstractVariable:
    """Base class for a variable that is the input to a function."""

    @property
    @abstractmethod
    def domain(self) -> AbstractDomain:
        """The domain that this variable is in."""

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check for equality."""

    @abstractmethod
    def __hash__(self) -> int:
        """Hash."""


class Variable(AbstractVariable):
    """A variable that is the input to a function."""

    _n = count(0)

    def __init__(self, domain: AbstractDomain, label: str | None = None):
        """Initialise."""
        if label is None:
            self._label = f"variable-{next(self._n)}"
        else:
            self._label - label
        self._domain = domain

    @property
    def label(self) -> str:
        """The label of this variable."""
        return self._label

    @property
    def domain(self) -> AbstractDomain:
        """The domain that this variable is in."""
        return self._domain

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, Variable) and other.label == self.label

    def __repr__(self) -> str:
        """Representation."""
        return f"Variable({self._label})"

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.Variable", self._label))


class FiniteElementVariable(AbstractVariable):
    """A variable that is the input to a function."""

    _n = count(0)

    def __init__(
        self, domain: AbstractFiniteElementDomain, label: str | None = None, reference: bool = False
    ):
        """Initialise."""
        if label is None:
            self._label = f"variable-{next(self._n)}"
        else:
            self._label = label
        self._domain = domain
        self._reference = reference

    @property
    def label(self) -> str:
        """The label of this variable."""
        return self._label

    @property
    def domain(self) -> AbstractDomain:
        """The domain that this variable is in."""
        return self._domain

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, FiniteElementVariable)
            and other.label == self.label
            and other.is_reference == self.is_reference
        )

    def __repr__(self) -> str:
        """Representation."""
        if self._reference:
            return f"FiniteElementVariable({self._label}, reference=True)"
        else:
            return f"FiniteElementVariable({self._label})"

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.FiniteElementVariable", self._label, self._reference))

    @property
    def is_reference(self) -> bool:
        """Check if this domain is on a reference cell."""
        return self._reference

    def to_reference(self) -> FiniteElementVariable:
        """Make a version of this variable on the reference cell."""
        return FiniteElementVariable(self._domain, self._label, True)

    def to_physical(self) -> FiniteElementVariable:
        """Make a version of this variable on physical cells."""
        return FiniteElementVariable(self._domain, self._label, True)


class AbstractFunction(AbstractExpression):
    """Base class for a function."""

    @property
    @abstractmethod
    def variable(self) -> AbstractVariable | None:
        """Get the variable that is this function's input."""

    @property
    def is_reference(self) -> bool:
        """Check if this function is on a reference cell."""
        if isinstance(self.variable, FiniteElementVariable):
            return self.variable.is_reference
        else:
            return False

    @abstractmethod
    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""

    @abstractmethod
    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""

    @property
    @abstractmethod
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""

    @property
    @abstractmethod
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        if self.is_reference:
            assert isinstance(self.function_space, AbstractReferenceMappedFunctionSpace)
            return self.function_space.elements[0].reference_value_shape
        else:
            return self.function_space.value_shape

    @property
    def is_cellwise_constant(self) -> bool:
        """Whether this function's value is the same everywhere on every cell.

        True when every element of the function space is known to lie in
        the degree-0 Lagrange space (lagrange_superdegree == 0) and,
        for a physical (non-reference) function, every element's
        reference map is known to preserve that constancy when pushed
        forward (see AbstractReferenceMap.preserves_constant_values).
        """
        assert isinstance(self.function_space, AbstractReferenceMappedFunctionSpace)
        elements = self.function_space.elements
        if any(e.lagrange_superdegree != 0 for e in elements):
            return False
        if self.is_reference:
            return True
        return all(e.reference_map.preserves_constant_values for e in elements)

    @property
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""
        return self.function_space.domain.cells[0].topological_dimension

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        if self.function_space.real_valued:
            return self
        else:
            return Re(self)

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        if self.function_space.real_valued:
            return zero(self.value_shape)
        else:
            return Im(self)


class Argument(AbstractFunction):
    """A function that is a dimension of the tensor to be assembled."""

    def __init__(
        self,
        space: AbstractFunctionSpace,
        component: int,
        is_reference: bool = False,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            space: The function space that this function lives in
            component: The component of the finite element tensor
                       to be assembled that this function represents
            is_reference: Is this argument's domain the reference cell?
            variable: The variable that is this argument's input
        """
        if variable is not None:
            if isinstance(variable, FiniteElementVariable):
                assert is_reference == variable.is_reference
            else:
                assert not is_reference
        self._space = space
        self._is_reference = is_reference
        self._variable = variable
        self._component = component

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._component, self._is_reference, variable)

    @property
    def is_reference(self) -> bool:
        """Check if this function is on a reference cell."""
        return self._is_reference

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @property
    def variable(self) -> AbstractVariable | None:
        """Get the variable that is this function's input."""
        return self._variable

    @property
    def component_index(self) -> int:
        """The component of the finite element tensor that this function represents."""
        return self._component

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self._component, self._is_reference, self._variable

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
                if isinstance(new, Argument) and self.variable is not None and new.variable is None:
                    return new.reconstruct_with_variable(self.variable)
                return new


class Coefficient(AbstractFunction):
    """A coefficient.

    A Coefficient represents a known function, such as a previous solution, a
    material property, or any other field supplied at assembly time.
    """

    _n = count(0)

    def __init__(
        self,
        space: AbstractFunctionSpace,
        coefficient_label: str | None = None,
        is_reference: bool = False,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            space: The function space that this function lives in
            coefficient_label: The label for this coefficient
            is_reference: Is this argument's domain the reference cell?
            variable: The variable that is this argument's input
        """
        if variable is not None:
            if isinstance(variable, FiniteElementVariable):
                assert is_reference == variable.is_reference
            else:
                assert not is_reference
        self._space = space
        self._is_reference = is_reference
        self._variable = variable
        if coefficient_label is None:
            self._label = f"coefficient-{next(self._n)}"
        else:
            self._label = coefficient_label

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._label, self._is_reference, variable)

    @property
    def is_reference(self) -> bool:
        """Check if this function is on a reference cell."""
        return self._is_reference

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @property
    def label(self) -> str:
        """The unique label of this coefficient."""
        return self._label

    @property
    def variable(self) -> AbstractVariable | None:
        """Get the variable that is this function's input."""
        return self._variable

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self._label, self._is_reference, self._variable

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        if not 0 <= index < self.domain_size:
            raise ValueError(
                f"Derivative index {index} out of range for domain size {self.domain_size}"
            )
        if self.is_cellwise_constant:
            # A cellwise constant's derivative is a plain zero
            # Tensor/RealScalar -- not itself an AbstractFunction -- so
            # this is the one place that distinction has to be cast away
            # rather than widening AbstractFunction.diff's own contract
            # (which would break chained .diff().diff() calls elsewhere,
            # eg test_basis_functions.py, whose intermediate values are
            # statically typed as bare AbstractFunction).
            return cast(AbstractFunction, zero(self.value_shape))
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
                and old.label == self.label
            ):
                if (
                    isinstance(new, Coefficient)
                    and self.variable is not None
                    and new.variable is None
                ):
                    return new.reconstruct_with_variable(self.variable)
                return new

    def pull_back_to_reference(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the reference cell."""
        if self.is_reference:
            raise ValueError("Cannot pull back function already defined on reference")
        assert isinstance(self._space, AbstractReferenceMappedFunctionSpace)
        if self._variable is None:
            return PushedForward(
                self._space.elements[0].reference_map,
                Coefficient(self._space, self._label, True, None),
            )
        else:
            assert isinstance(self._variable, FiniteElementVariable)
            return PushedForward(
                self._space.elements[0].reference_map,
                Coefficient(self._space, self._label, True, self._variable.to_reference()),
            )


class TestFunction(Argument):
    """A test function."""

    __test__ = False

    def __init__(
        self,
        space: AbstractFunctionSpace,
        is_reference: bool = False,
        variable: AbstractVariable | None = None,
    ):
        """Initialise."""
        super().__init__(space, 0, is_reference, variable)

    def __repr__(self) -> str:
        """Representation."""
        return f"TestFunction({self.variable!r})"

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._is_reference, variable)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self.is_reference, self.variable

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()

    def pull_back_to_reference(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the reference cell."""
        if self.is_reference:
            raise ValueError("Cannot pull back function already defined on reference")
        assert isinstance(self._space, AbstractReferenceMappedFunctionSpace)
        if self.variable is None:
            return PushedForward(
                self._space.elements[0].reference_map,
                TestFunction(self._space, True, None),
            )
        else:
            assert isinstance(self.variable, FiniteElementVariable)
            return PushedForward(
                self._space.elements[0].reference_map,
                TestFunction(self._space, True, self.variable.to_reference()),
            )


class TrialFunction(Argument):
    """A trial function."""

    def __init__(
        self,
        space: AbstractFunctionSpace,
        is_reference: bool = False,
        variable: AbstractVariable | None = None,
    ):
        """Initialise."""
        super().__init__(space, 1, is_reference, variable)

    def __repr__(self) -> str:
        """Representation."""
        return f"TrialFunction({self.variable!r})"

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._is_reference, variable)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self.is_reference, self.variable

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()

    def pull_back_to_reference(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the reference cell."""
        if self.is_reference:
            raise ValueError("Cannot pull back function already defined on reference")
        assert isinstance(self._space, AbstractReferenceMappedFunctionSpace)
        if self.variable is None:
            return PushedForward(
                self._space.elements[0].reference_map,
                TrialFunction(self._space, True, None),
            )
        else:
            assert isinstance(self.variable, FiniteElementVariable)
            return PushedForward(
                self._space.elements[0].reference_map,
                TrialFunction(self._space, True, self.variable.to_reference()),
            )


def create_variable(domain: AbstractDomain) -> AbstractVariable:
    """Create a new variable in a domain."""
    if isinstance(domain, AbstractFiniteElementDomain):
        return FiniteElementVariable(domain)
    else:
        return Variable(domain)
