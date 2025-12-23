# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Finite element function spaces."""

from abc import ABC, abstractmethod

from uflx.domain import AbstractDomain
from uflx.finite_element import AbstractFiniteElement


class AbstractFunctionSpace(ABC):
    """Abstract base class for a function space."""

    @property
    @abstractmethod
    def domains(self) -> list[AbstractDomain]:
        """Domains included in this function space."""


class FunctionSpace(AbstractFunctionSpace):
    """Function space."""

    def __init__(self, domains_and_elements: list[tuple[AbstractDomain, AbstractFiniteElement]]):
        """Initialise."""
        self._d_and_e = domains_and_elements

    @property
    def domains_and_elements(self) -> list[tuple[AbstractDomain, AbstractFiniteElement]]:
        """Domains and elements included in this function space."""
        return self._d_and_e

    @property
    def domains(self) -> list[AbstractDomain]:
        """Domains included in this function space."""
        return [i[0] for i in self.domains_and_elements]

    @property
    def elements(self) -> list[AbstractFiniteElement]:
        """Elements included in this function space."""
        return [i[1] for i in self.domains_and_elements]
