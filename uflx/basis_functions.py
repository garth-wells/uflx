# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Finite element basis functions."""

from __future__ import annotations

from abc import abstractmethod
from typing import Any

from uflx.expressions import AbstractExpression, Im, Re
from uflx.finite_elements import AbstractFiniteElement, AbstractReferenceMappedFiniteElement
from uflx.function_spaces import AbstractFunctionSpace, AbstractReferenceMappedFunctionSpace
from uflx.functions import AbstractFunction
from uflx.graphs import GraphNode
from uflx.points import AbstractPoint
from uflx.tensors import zero
from uflx.utils import flatten


class AbstractEvaluatedBasisFunction(AbstractFunction):
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

    @abstractmethod
    def diff(self, index: int) -> AbstractEvaluatedBasisFunction:
        """Take a derivative of this function."""

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        if self.element.real_valued:
            return self
        else:
            return Re(self)

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        if self.element.real_valued:
            return zero(self.value_shape)
        else:
            return Im(self)


class EvaluatedBasisFunction(AbstractEvaluatedBasisFunction):
    """A basis function evaluated at a point."""

    def __init__(
        self,
        space: AbstractReferenceMappedFunctionSpace,
        basis_index: int | str,
        point: AbstractPoint,
        is_reference: bool,
        element_index: int | None = None,
        derivative: tuple[int, ...] | None = None,
        component: int | None = None,
    ):
        """Initialise."""
        self._space = space
        if element_index is None:
            if len(space.elements) > 1:
                raise ValueError(
                    "Basis functions in spaces with more than one element "
                    "must be given an element index."
                )
            self._element = space.elements[0]
        else:
            self._element = space.elements[element_index]
        self._element_index = element_index
        self._basis_index = basis_index
        self._point = point
        if derivative is None:
            self._derivative = tuple(0 for _ in range(self._element.cell.topological_dimension))
        else:
            self._derivative = derivative
        if component is None and self._element.reference_value_size == 1:
            self._component: int | None = 0
        else:
            self._component = component
        self._is_reference = is_reference

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @property
    def is_reference(self) -> bool:
        """Is this function's domain the reference cell?"""
        return self._is_reference

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
        return self._point.index

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        if self._component is None:
            if self.is_reference:
                assert isinstance(self.element, AbstractReferenceMappedFiniteElement)
                return self.element.reference_value_shape
            else:
                return self.element.physical_value_shape(self._point.dim)
        else:
            return ()

    def __repr__(self):
        """Representation."""
        repr = (
            "EvaluatedBasisFunction("
            f"{self._space!r}, {self._basis_index}, {self._point!r}, {self.is_reference}"
        )
        if self._element_index is not None:
            repr += f", {self._element_index}"
        if self._derivative is not None:
            repr += f", derivative={self._derivative}"
        if self._component is not None:
            repr += f", component={self._component}"
        repr += ")"
        return repr

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (
            self._space,
            self._basis_index,
            self._point,
            self.is_reference,
            self._element_index,
            self._derivative,
            self._component,
        )

    @property
    def derivative(self) -> tuple[int, ...]:
        """The number of derivatives in each coordinate direction."""
        return self._derivative

    @property
    def component_index(self) -> int | None:
        """The (flattened) component index of the basis function."""
        return self._component

    def diff(self, index: int) -> EvaluatedBasisFunction:
        """Take a derivative of this function."""
        return EvaluatedBasisFunction(
            self._space,
            self._basis_index,
            self._point,
            self.is_reference,
            self._element_index,
            tuple(d + 1 if i == index else d for i, d in enumerate(self._derivative)),
            self._component,
        )

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return EvaluatedBasisFunction(
            self._space,
            self._basis_index,
            self._point,
            self.is_reference,
            self._element_index,
            self._derivative,
            flatten(indices, self.value_shape),
        )

    @property
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""
        return self.element.cell.topological_dimension
