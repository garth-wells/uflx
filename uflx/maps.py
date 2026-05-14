"""Push forward and pull back maps.

These maps are uses to map function values between reference cells and physical cells
"""

from abc import ABC, abstractmethod


class AbstractReferenceMap(ABC):
    """Abstract base class for reference maps."""

    @abstractmethod
    def push_forward_symbolic(self, function):
        """Map function values from a reference cell to a physical cell."""

    @abstractmethod
    def pull_back_symbolic(self, function):
        """Map function values from a physical cell to a reference cell."""


class IdentityReferenceMap(AbstractReferenceMap):
    """Indentity map."""

    def push_forward_symbolic(self, function):
        """Map function values from a reference cell to a physical cell."""
        return function

    def pull_back_symbolic(self, function):
        """Map function values from a physical cell to a reference cell."""
        return function
