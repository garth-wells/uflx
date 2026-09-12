"""Graph Nodes representing code structures."""

from typing import Any

from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode

from uflx_codegeneration import symbols
from uflx_codegeneration.c import GenerateC
from uflx_codegeneration.utils import indented


def flatten_component(
    indices: tuple[int | str, ...],
    shape: tuple[int, ...],
    bracketed: bool = False,
) -> str:
    """Flatten the component in an array access."""
    assert len(indices) == len(shape)
    if len(indices) == 1:
        return str(indices[0])

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
        shape: tuple[int, ...],
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

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


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

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class Declare:
    """Declare and initialise a scalar variable."""

    def __init__(self, dtype: str, variable: str, value: float | int):
        """Initalise."""
        self.dtype = dtype
        self.variable = variable
        self.value = value

    def __repr__(self):
        """Representation."""
        return f"Declare({self.dtype}, {self.variable}, {self.value})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.dtype, self.variable, self.value

    def generate_c(self) -> str:
        """Generate code for this object."""
        return f"{self.dtype} {self.variable} = {self.value};"


class AccumulateToVariable:
    """Add the value of an expression into an existing scalar variable."""

    def __init__(self, variable: str, body: GraphNode):
        """Initalise."""
        self.variable = variable
        self.body = body

    def __repr__(self):
        """Representation."""
        return f"AccumulateToVariable({self.variable})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.variable, self.body

    def generate_c(self) -> str:
        """Generate code for this object."""
        assert isinstance(self.body, GenerateC)
        return f"{self.variable} += {self.body.generate_c()};"


class Return:
    """Return a value from a function."""

    def __init__(self, body: GraphNode):
        """Initalise."""
        self.body = body

    def __repr__(self):
        """Representation."""
        return f"Return({self.body})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.body,)

    def generate_c(self) -> str:
        """Generate code for this object."""
        assert isinstance(self.body, GenerateC)
        return f"return {self.body.generate_c()};"


class Block:
    """A sequence of statements, executed in order.

    This is standard (C99) C's equivalent of a compound statement: unlike a
    GNU/Clang statement-expression, a Block cannot itself be used as a value
    inside a larger expression, only as the body (or part of the body) of a
    function or loop.
    """

    def __init__(self, statements: tuple[Any, ...]):
        """Initalise."""
        self.statements = statements

    def __repr__(self):
        """Representation."""
        return f"Block({self.statements!r})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {s for s in self.statements if isinstance(s, GraphNode)}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self.statements,)

    def generate_c(self) -> str:
        """Generate code for this object."""
        parts = []
        for statement in self.statements:
            assert isinstance(statement, GenerateC)
            parts.append(statement.generate_c())
        return "\n".join(parts)


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

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")
