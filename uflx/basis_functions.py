# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Finite element basis functions."""

from __future__ import annotations

from abc import abstractmethod
from typing import Any

from uflx.expressions import AbstractExpression
from uflx.finite_elements import AbstractFiniteElement, AbstractReferenceMappedFiniteElement
from uflx.function_spaces import AbstractFunctionSpace
from uflx.functions import AbstractPhysicalFunction, AbstractReferenceFunction
from uflx.graphs import GraphNode
from uflx.points import AbstractPoint, AbstractPointInSet
from uflx.utils import flatten


class AbstractEvaluatedReferenceBasisFunction(AbstractReferenceFunction):
    """Base class for a basis function evaluated at a point on the reference cell."""

    @property
    @abstractmethod
    def element(self) -> AbstractFiniteElement:
        """The finite element containing this basis function."""

    @property
    @abstractmethod
    def basis_index(self) -> int | str:
        """The index of the basis function."""

    @property
    @abstractmethod
    def point_index(self) -> int | str:
        """The index of the point in the set of points."""

    @property
    @abstractmethod
    def point(self) -> AbstractPoint:
        """The point at which the function is evaluated."""

    @property
    @abstractmethod
    def derivative(self) -> tuple[int, ...]:
        """The number of derivatives in each coordinate direction."""

    @property
    @abstractmethod
    def component_index(self) -> int | None:
        """The (flattened) component index of the basis function."""


class EvaluatedReferenceBasisFunction(AbstractEvaluatedReferenceBasisFunction):
    """A basis function evaluated at a point on the reference cell."""

    def __init__(
        self,
        element: AbstractFiniteElement,
        basis_index: int | str,
        point: AbstractPoint,
        derivative: tuple[int, ...] | None = None,
        component: int | None = None,
    ):
        """Initialise."""
        self._element = element
        self._basis_index = basis_index
        self._point = point
        if derivative is None:
            self._derivative = tuple(0 for _ in range(element.cell.topological_dimension))
        else:
            self._derivative = derivative
        if (
            component is None
            and isinstance(element, AbstractReferenceMappedFiniteElement)
            and element.reference_value_size == 1
        ):
            self._component: int | None = 0
        else:
            self._component = component

    @property
    def point(self) -> AbstractPoint:
        """The point at which the function is evaluated."""
        return self._point

    @property
    def element(self) -> AbstractFiniteElement:
        """The finite element containing this basis function."""
        return self._element

    @property
    def basis_index(self) -> int | str:
        """The index of the basis function."""
        return self._basis_index

    @property
    def point_index(self) -> int | str:
        """The index of the point in the set of points."""
        if isinstance(self._point, AbstractPointInSet):
            return self._point.index
        return 0

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def __repr__(self):
        """Representation."""
        repr = (
            "EvaluatedReferenceBasisFunction("
            f"{self._element!r}, {self._basis_index}, {self._point!r}"
        )
        if self._derivative is not None:
            repr += f", {self._derivative}"
        if self._component is not None:
            repr += f", {self._component}"
        repr += ")"
        return repr

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._element, self._basis_index, self._point, self._derivative, self._component

    @property
    def derivative(self) -> tuple[int, ...]:
        """The number of derivatives in each coordinate direction."""
        return self._derivative

    @property
    def component_index(self) -> int | None:
        """The (flattened) component index of the basis function."""
        return self._component

    def diff(self, index: int) -> AbstractReferenceFunction:
        """Take a derivative of this function."""
        return EvaluatedReferenceBasisFunction(
            self._element,
            self._basis_index,
            self._point,
            tuple(d + 1 if i == index else d for i, d in enumerate(self._derivative)),
            self._component,
        )

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return EvaluatedReferenceBasisFunction(
            self._element,
            self._basis_index,
            self._point,
            self._derivative,
            flatten(indices, self.value_shape),
        )

    @property
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""
        return self.element.cell.topological_dimension


class AbstractEvaluatedPhysicalBasisFunction(AbstractPhysicalFunction):
    """Base class for a basis function evaluated at a point on a physical cell."""

    @property
    @abstractmethod
    def element(self) -> AbstractFiniteElement:
        """The finite element containing this basis function."""

    @property
    @abstractmethod
    def basis_index(self) -> int | str:
        """The index of the basis function."""

    @property
    @abstractmethod
    def point_index(self) -> int | str:
        """The index of the point in the set of points."""

    @property
    @abstractmethod
    def point(self) -> AbstractPoint:
        """The point at which the function is evaluated."""


class EvaluatedPhysicalBasisFunction(AbstractEvaluatedPhysicalBasisFunction):
    """A basis function evaluated at a point on the physical cell."""

    def __init__(
        self,
        function_space: AbstractFunctionSpace,
        element: AbstractFiniteElement,
        basis_index: int | str,
        point: AbstractPoint,
        derivative: tuple[int, ...] | None = None,
        component: int | None = None,
    ):
        """Initalise."""
        self._function_space = function_space
        self._element = element
        self._basis_index = basis_index
        self._point = point
        if derivative is None:
            self._derivative = tuple(0 for _ in range(element.cell.topological_dimension))
        else:
            self._derivative = derivative
        if (
            component is None
            and isinstance(element, AbstractReferenceMappedFiniteElement)
            and element.reference_value_size == 1
        ):
            self._component: int | None = 0
        else:
            self._component = component

    @property
    def point(self) -> AbstractPoint:
        """The point at which the function is evaluated."""
        return self._point

    @property
    def element(self) -> AbstractFiniteElement:
        """The finite element containing this basis function."""
        return self._element

    @property
    def basis_index(self) -> int | str:
        """The index of the basis function."""
        return self._basis_index

    @property
    def point_index(self) -> int | str:
        """The index of the point in the set of points."""
        if isinstance(self._point, AbstractPointInSet):
            return self._point.index
        return 0

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def __repr__(self):
        """Representation."""
        repr = (
            "EvaluatedPhysicalBasisFunction("
            f"{self._element!r}, {self._basis_index}, {self._point!r}"
        )
        if self._derivative is not None:
            repr += f", {self._derivative}"
        if self._component is not None:
            repr += f", {self._component}"
        repr += ")"
        return repr

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._element, self._basis_index, self._point, self._derivative, self._component

    @property
    def derivative(self) -> tuple[int, ...]:
        """The number of derivatives in each coordinate direction."""
        return self._derivative

    @property
    def component_index(self) -> int | None:
        """The (flattened) component index of the basis function."""
        if self._component is None:
            raise ValueError("EvaluatedBasisFunction is not an evaluation of a single component")
        return self._component

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._function_space

    @property
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""
        return self.element.cell.topological_dimension

    def diff(self, index: int) -> AbstractReferenceFunction:
        """Take a derivative of this function."""
        return EvaluatedReferenceBasisFunction(
            self._element,
            self._basis_index,
            self._point,
            tuple(d + 1 if i == index else d for i, d in enumerate(self._derivative)),
            self._component,
        )

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return EvaluatedReferenceBasisFunction(
            self._element,
            self._basis_index,
            self._point,
            self._derivative,
            flatten(indices, self.value_shape),
        )
