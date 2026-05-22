"""Graph Nodes representing code structures."""

from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode
from uflx.scalars import AbstractInteger

from uflx_codegeneration import symbols
from uflx_codegeneration.c import GenerateC
from uflx_codegeneration.utils import indented


def flatten_component(
    indices: tuple[int | str, ...],
    shape: tuple[int | AbstractInteger, ...],
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

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.variable, self.start, self.end, self.body

    def generate_c(self) -> str:
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
        shape: tuple[int | AbstractInteger, ...],
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

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.component, self.shape, self.body

    def __repr__(self):
        """Representation."""
        return f"AddToLocalTensor({self.component})"

    def generate_c(self) -> str:
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

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.array, self.index

    def __repr__(self):
        """Representation."""
        return f"{self.array}[{','.join(str(i) for i in self.index)}]"

    def generate_c(self) -> str:
        """Generate code for this object."""
        return f"{self.array}[" + "][".join(f"{i}" for i in self.index) + "]"


class FunctionCall(AbstractExpression):
    """A call to a function."""

    def __init__(self, function: str, *inputs: Any):
        """Initalise."""
        self.function = function
        self.inputs = inputs

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {i for i in self.inputs if isinstance(i, GraphNode)}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.function, *self.inputs

    def __repr__(self):
        """Representation."""
        return f"FunctionCall({self.function}, (" + ", ".join(f"{i!r}" for i in self.inputs) + "))"

    def generate_c(self) -> str:
        """Generate code for this object."""
        return (
            f"{self.function}("
            + ", ".join(i.generate_c() if isinstance(i, GenerateC) else f"{i}" for i in self.inputs)
            + ")"
        )


class Variable(AbstractExpression):
    """A variable."""

    def __init__(self, dtype: str, variable: str):
        """Initialise."""
        self._dtype = dtype
        self._variable = variable

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._dtype, self._variable

    def __repr__(self):
        """Representation."""
        return f"Variable({self._dtype}, {self._variable})"

    def generate_c(self) -> str:
        """Generate code for this object."""
        return self._variable
