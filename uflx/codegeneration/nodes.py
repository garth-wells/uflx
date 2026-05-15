"""Graph Nodes representing code structures."""

from typing import Self

from uflx.codegeneration import symbols
from uflx.codegeneration.c import GenerateC
from uflx.expressions import AbstractExpression
from uflx.finite_elements import Dimension
from uflx.graphs import GraphNode
from uflx.utils import indented


def flatten_component(
    indices: tuple[int | str, ...],
    shape: tuple[int | Dimension, ...],
    bracketed: bool = False,
):
    """Flatten the component in an array access."""
    assert len(indices) == len(shape)
    if len(indices) == 1:
        return indices[0]

    component = (
        flatten_component(indices[:-1], shape[:-1], True) + f" * {shape[-1]} + {indices[-1]}"
    )
    if bracketed:
        return f"({component})"
    else:
        return f"{component}"


class Loop:
    """A for loop."""

    def __init__(self, variable: str, start: int | str, end: int | str, body: GraphNode):
        """Initalise."""
        self.variable = variable
        self.start = start
        self.end = end
        self.body = body

    def __repr__(self):
        """Representation."""
        return f"Loop({self.variable}, {self.start}, {self.end})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self.__class__(
            self.variable, self.start, self.end, replacements.get(self.body, self.body)
        )

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.body, GenerateC)
        return (
            f"for (int {self.variable}={self.start}; {self.variable}!={self.end}; "
            f"++{self.variable})\n"
            "{\n" + indented(self.body.generate_c(), 2) + "\n}"
        )


class AddToLocalTensor:
    """Add to an entry in the local tensor for the current cell."""

    def __init__(
        self,
        component: tuple[int | str, ...],
        shape: tuple[int | Dimension, ...],
        body: AbstractExpression,
    ):
        """Initalise."""
        self.component = component
        self.shape = shape
        self.body = body

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        body = replacements.get(self.body, self.body)
        assert isinstance(body, AbstractExpression)
        return self.__class__(self.component, self.shape, body)

    def __repr__(self):
        """Representation."""
        return f"AddToLocalTensor({self.component})"

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        assert isinstance(self.body, GenerateC)
        return (
            f"{symbols.local_tensor}["
            + flatten_component(self.component, self.shape)
            + "] += "
            + self.body.generate_c()
            + ";"
        )


class ArrayEntry(AbstractExpression):
    """A single item in an array."""

    def __init__(self, array: str, index: tuple[int | str, ...]):
        """Initalise."""
        self.array = array
        self.index = index

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self.__class__(self.array, self.index)

    def __repr__(self):
        """Representation."""
        return f"ArrayEntry({self.array}, {self.index})"

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        return f"{self.array}[" + "][".join(f"{i}" for i in self.index) + "]"
