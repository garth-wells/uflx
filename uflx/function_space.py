"""Function space."""

from abc import ABC, abstractmethod

from uflx.finite_element import AbstractFiniteElement
from uflx.domain import AbstractDomain


class AbstractFunctionSpace(ABC):
    """Abstract base class for a function space."""

    @property
    @abstractmethod
    def domains(self) -> list[AbstractDomain]:
        """Domains included in this function space."""

    @property
    @abstractmethod
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the function space."""


class FunctionSpace(AbstractFunctionSpace):
    """Function space."""

    def __init__(self, domains_and_elements: list[tuple[AbstractDomain, AbstractFiniteElement]]):
        """Initialise."""
        self._d_and_e = domains_and_elements

        gdim = domains_and_elements[0][0].geometric_dimension
        shape = domains_and_elements[0][1].physical_value_shape(gdim)
        for domain, element in domains_and_elements[1:]:
            if domain.geometric_dimension != gdim:
                raise ValueError(
                    "Domains in a functions space must have the same geometric dimension."
                )
            if element.physical_value_shape(gdim) != shape:
                raise ValueError(
                    "Elements in a functions space must have the same physical value shape."
                )

    @property
    def domains_and_elements(self) -> list[tuple[AbstractDomain, AbstractFiniteElement]]:
        return self._d_and_e

    @property
    def domains(self) -> list[AbstractDomain]:
        """Domains included in this function space."""
        return [i[0] for i in self.domains_and_elements]

    @property
    def elements(self) -> list[AbstractFiniteElement]:
        """Elements included in this function space."""
        return [i[1] for i in self.domains_and_elements]

    @property
    def value_shape(self) -> tuple[int, ...]:
        return self.elements[0].physical_value_shape(self.domains[0].geometric_dimension)
