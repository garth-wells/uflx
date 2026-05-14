"""Code generation."""

from typing import Protocol, Self

import networkx as nx
import numpy as np

from uflx.domains import AbstractCoordinateElement
from uflx.function_spaces import AbstractReferenceMappedFunctionSpace
from uflx.functions import Argument
from uflx.graphs import Graph, GraphNode, RepresentedByGraph, generate_graph, is_dag
from uflx.graphs.algorithms import replace
from uflx.integrals import AbstractIntegral, dx
from uflx.operators import Abs, Add, Conj, Mult, Subtract

# TODO: import these from codegeneration.symbols
local_tensor = "A"
coordinate_dofs = "coordinate_dofs"


class ConvertToCCode(Protocol):
    """Protocol for Objects that can be converted to C code."""

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""


def reconstruct(object, args, replacements):
    """Reconstruct an object from its arguements with replacements made."""
    if all(a not in replacements for a in args):
        return object
    return object.__class__(*(replacements.get(a, a) for a in args))


def print_node(graph: Graph, node: GraphNode, indentation: int = 0):
    """Print a graph using the node as the root node."""
    print(" " * (2 * indentation) + f"{node!r}")
    for next in graph.successors(node):
        print_node(graph, next, indentation + 1)


def print_graph(graph: Graph):
    """Print a graph."""
    print_node(graph, graph.root)


class QuadratureRule:
    """A quadrature rule."""

    def __init__(self, points, weights):
        """Initialise."""
        self.points = points
        self.weights = weights

    @property
    def npoints(self):
        """The number of points in the quadrature rule."""
        return len(self.weights)


class QuadraturePoint:
    """A point in a quadrature rule."""

    def __init__(self, rule, index):
        """Initalise."""
        self.rule = rule
        self.index = index

    def __repr__(self):
        """Representation."""
        return f"QuadraturePoint({self.index})"


class QuadratureWeight:
    """A weight in a quadrature rule."""

    def __init__(self, rule, index):
        """Initalise."""
        self.rule = rule
        self.index = index

    def __repr__(self):
        """Representation."""
        return f"QuadratureWeight({self.index})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self


class JacobianDeterminant:
    """The determinant of the Jacobian."""

    def __init__(self, domain, point):
        """Initalise."""
        self.domain = domain
        self.point = point

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self


class QuadratureLoop:
    """A loop over the points in a quadrature rule."""

    def __init__(self, body, rule, variable):
        """Initalise."""
        self.body = body
        self.rule = rule
        self.variable = variable

    def __repr__(self):
        """Representation."""
        return f"QuadratureLoop({self.variable})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return reconstruct(self, (self.body, self.rule, self.variable), replacements)

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        return (
            f"for (int {self.variable}=0; {self.variable}!={self.rule.npoints}; "
            f"++{self.variable})\n"
            "{\n" + indented(self.body.generate_c(), 2) + "\n}"
        )


class Scalar:
    """A scalar."""

    def __init__(self, value):
        """Initalise."""
        self.value = value

    def __repr__(self):
        """Representation."""
        return f"{self.value}"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        return f"{self.value}"


def indented(code, spaces):
    """Add indentation to a block of code."""
    return "\n".join(" " * spaces + line for line in code.split("\n"))


class Loop:
    """A for loop."""

    def __init__(self, variable: str, start: int, end: int, body):
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
        return reconstruct(self, (self.variable, self.start, self.end, self.body), replacements)

    def generate_c(self, bracketed: bool = False) -> list[str]:
        """Generate code for this object."""
        return (
            f"for (int {self.variable}={self.start}; {self.variable}!={self.end}; "
            f"++{self.variable})\n"
            "{\n" + indented(self.body.generate_c(), 2) + "\n}"
        )


class EvaluatedBasisFunction:
    """A basis function evaluated at a point."""

    def __init__(self, element, basis_index, point):
        """Initalise."""
        self.element = element
        self.basis_index = basis_index
        self.point = point

    def __repr__(self):
        """Representation."""
        return f"EvaluatedBasisFunction({self.element!r}, {self.basis_index}, {self.point!r})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self


class EvaluatedBasisFunctionDerivative:
    """A derivative of a basis function evaluated at a point."""

    def __init__(self, element, basis_index, point, derivative):
        """Initalise."""
        self.element = element
        self.basis_index = basis_index
        self.point = point
        self.derivative = derivative

    def __repr__(self):
        """Representation."""
        return (
            f"EvaluatedBasisFunctionDerivative({self.element!r}, {self.basis_index}, "
            f"{self.point!r}, {self.derivative})"
        )

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self


class PushedForward:
    """A function on a reference cell that has been mapped to a physical cell."""

    def __init__(self, map, function):
        """Initalise."""
        self.map = map
        self.function = function

    def __repr__(self):
        """Representation."""
        return f"PushedForward({self.map})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.function}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return reconstruct(self, (self.map, self.function), replacements)


class VariableNamer:
    """Class for naming variables to ensure that temporary variables do not have clashing names."""

    def __init__(self):
        """Initalise."""
        self.i = -1
        self.n = -1
        self.fe_i = -1
        self.qr_i = -1

    def variable(self):
        """Get a new variable name."""
        chars = ["i", "j", "k"]
        self.i += 1
        if self.i == len(chars):
            self.i = 0
            self.n += 1
        if self.n == -1:
            return chars[self.i]
        else:
            return f"{chars[self.i]}{self.n}"

    def finite_element_table(self):
        """Get a new finite element table name."""
        self.fe_i += 1
        return f"FE{self.fe_i}"

    def quadrature_table(self):
        """Get a new quadrature weight table name."""
        self.qr_i += 1
        return f"QW{self.qr_i}"


global_variable_namer = VariableNamer()


def flatten_component(indices, shape, bracketed=False):
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


class AddToLocalTensor:
    """Add to an entry in the local tensor for the current cell."""

    def __init__(self, component, shape, body):
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
        return reconstruct(self, (self.component, self.shape, self.body), replacements)

    def __repr__(self):
        """Representation."""
        return f"AddToLocalTensor({self.component})"

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        return (
            f"{local_tensor}["
            + flatten_component(self.component, self.shape)
            + "] += "
            + self.body.generate_c()
            + ";"
        )


class ArrayEntry:
    """A single item in an array."""

    def __init__(self, array, index):
        """Initalise."""
        self.array = array
        self.index = index

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self

    def __repr__(self):
        """Representation."""
        return f"ArrayEntry({self.array}, {self.index})"

    def generate_c(self, bracketed: bool = False) -> str:
        """Generate code for this object."""
        return f"{self.array}[" + "][".join(f"{i}" for i in self.index) + "]"


def apply_push_forwards(
    graph: Graph,
) -> Graph:
    """Apply push forward maps to functions."""
    return replace(
        graph,
        {
            node: node.map.push_forward_symbolic(node.function)
            for node in graph
            if isinstance(node, PushedForward)
        },
    )


def convert_complex_to_real(
    graph: Graph,
) -> Graph:
    """Take the real part of all complex values."""
    return replace(graph, {node: node.init_arg for node in graph if isinstance(node, Conj)})


def integrals_to_quadrature(
    graph: Graph,
    rules: dict[RepresentedByGraph, QuadratureRule],
    variable_namer=global_variable_namer,
) -> Graph:
    """Replace integrals with quadrature."""
    updated_nodes = {}
    to_replace = {}

    for node in reversed(list(nx.topological_sort(graph))):
        if isinstance(node, AbstractIntegral):
            tensor_shape_components = {}

            if node.measure._codim != 0 or node.measure._boundary_only:
                raise NotImplementedError("Only codim 0 integrals supported for now")

            arguments = []
            for i in nx.descendants(graph, node):
                if isinstance(i, Argument):
                    arguments.append(i)
            for i in arguments:
                i_space = i.function_space
                if not isinstance(i_space, AbstractReferenceMappedFunctionSpace):
                    raise NotImplementedError(
                        "Code generation only implemented for reference mapped spaces"
                    )
                if len(i_space.elements) != 1:
                    raise NotImplementedError(
                        "Code generation currently only implemented for spaces with "
                        "exactly one element"
                    )
                tensor_shape_components[i.component] = i_space.elements[0].dim
            tensor_shape = tuple(
                tensor_shape_components.get(i, 1)
                for i in range(
                    min(tensor_shape_components.keys()), max(tensor_shape_components.keys()) + 1
                )
            )
            variables = tuple(
                0 if component == 1 else variable_namer.variable() for component in tensor_shape
            )

            qvariable = variable_namer.variable()

            for a in arguments:
                point = QuadraturePoint(rules[node.measure], variables[a.component])
                to_replace[a] = PushedForward(
                    a.function_space.elements[0].map_type,
                    EvaluatedBasisFunction(a.function_space.elements[0], qvariable, point),
                )

            domain = arguments[0].function_space.domain
            for a in arguments:
                assert a.function_space.domain == domain

            integrand = Mult(
                Mult(
                    QuadratureWeight(rules[node.measure], variables[a.component]),
                    Abs(
                        JacobianDeterminant(domain, QuadraturePoint(rules[node.measure], qvariable))
                    ),
                ),
                node.integrand,
            )

            body = AddToLocalTensor(variables, tensor_shape, integrand)

            qloop = QuadratureLoop(body, rules[node.measure], qvariable)

            next = qloop
            for variable, count in zip(variables[::-1], tensor_shape[::-1]):
                loop = Loop(variable, 0, count, next)
                next = loop

            updated_nodes[node] = next

    new_graph = generate_graph(updated_nodes.get(graph.root, graph.root))
    new_graph = replace(new_graph, to_replace)
    return new_graph


def tabulate_finite_elements(
    graph,
    variable_namer=global_variable_namer,
):
    """Generate tables of values for finite elements that need to be evaluated."""
    table_map = {}
    tables = {}
    to_replace = {}
    for node in graph:
        if isinstance(node, EvaluatedBasisFunction):
            if not isinstance(node.point, QuadraturePoint):
                raise NotImplementedError()

            id = (node.element, node.point.rule, 0)
            if id not in table_map:
                name = variable_namer.finite_element_table()
                table_map[id] = name
                tables[name] = node.element.tabulate(
                    node.point.rule.points,
                    tuple(0 for _ in range(node.element.cell.topological_dimension)),
                )
            to_replace[node] = ArrayEntry(table_map[id], (node.point.index, node.basis_index))

        if isinstance(node, EvaluatedBasisFunctionDerivative):
            if not isinstance(node.point, QuadraturePoint):
                raise NotImplementedError()

            id = (node.element, node.point.rule, node.derivative)
            if id not in table_map:
                name = variable_namer.finite_element_table()
                table_map[id] = name
                tables[name] = node.element.tabulate(node.point.rule.points, node.derivative)
            to_replace[node] = ArrayEntry(table_map[id], (node.point.index, node.basis_index))

    return tables, replace(graph, to_replace)


def tabulate_quadrature(
    graph,
    variable_namer=global_variable_namer,
):
    """Generate tables of values for quadrature rules."""
    table_map = {}
    tables = {}
    to_replace = {}
    for node in graph:
        if isinstance(node, QuadratureWeight):
            if node.rule not in table_map:
                name = variable_namer.quadrature_table()
                table_map[node.rule] = name
                tables[name] = np.array(node.rule.weights)
            to_replace[node] = ArrayEntry(table_map[node.rule], (node.index,))

    return tables, replace(graph, to_replace)


def expand_jacobians(
    graph,
    variable_namer=global_variable_namer,
) -> Graph:
    """Replace jacobians with evaluations of the derivatives of finite elements."""
    to_replace = {}

    for node in graph:
        if isinstance(node, JacobianDeterminant):
            if not isinstance(node.domain, AbstractCoordinateElement):
                raise NotImplementedError()
            if len(node.domain.elements) > 1:
                raise NotImplementedError()
            (element,) = node.domain.elements
            tdim = element.cell.topological_dimension
            gdim = node.domain.geometric_dimension

            if tdim == 0 and gdim == 0:
                to_replace[node] = Scalar(1.0)
            elif tdim == 2 and gdim == 2:
                j00 = Mult(
                    ArrayEntry(coordinate_dofs, (0,)),
                    EvaluatedBasisFunctionDerivative(element, 0, node.point, (1, 0)),
                )
                j01 = Mult(
                    ArrayEntry(coordinate_dofs, (0,)),
                    EvaluatedBasisFunctionDerivative(element, 0, node.point, (0, 1)),
                )
                j10 = Mult(
                    ArrayEntry(coordinate_dofs, (1,)),
                    EvaluatedBasisFunctionDerivative(element, 0, node.point, (1, 0)),
                )
                j11 = Mult(
                    ArrayEntry(coordinate_dofs, (1,)),
                    EvaluatedBasisFunctionDerivative(element, 0, node.point, (0, 1)),
                )

                for i in range(1, element.dim):
                    j00 = Add(
                        j00,
                        Mult(
                            ArrayEntry(coordinate_dofs, (tdim * i,)),
                            EvaluatedBasisFunctionDerivative(element, i, node.point, (1, 0)),
                        ),
                    )
                    j01 = Add(
                        j01,
                        Mult(
                            ArrayEntry(coordinate_dofs, (tdim * i,)),
                            EvaluatedBasisFunctionDerivative(element, i, node.point, (0, 1)),
                        ),
                    )
                    j10 = Add(
                        j10,
                        Mult(
                            ArrayEntry(coordinate_dofs, (tdim * i + 1,)),
                            EvaluatedBasisFunctionDerivative(element, i, node.point, (1, 0)),
                        ),
                    )
                    j11 = Add(
                        j11,
                        Mult(
                            ArrayEntry(coordinate_dofs, (tdim * i + 1,)),
                            EvaluatedBasisFunctionDerivative(element, i, node.point, (0, 1)),
                        ),
                    )

                to_replace[node] = Subtract(Mult(j00, j11), Mult(j01, j10))
            else:
                raise NotImplementedError()
    return replace(graph, to_replace)


def c_table(table: np.ndarray) -> str:
    """Convert a numpy array to C."""
    if len(table.shape) == 1:
        return "{" + ", ".join(f"{i}" for i in table) + "}"
    return "{" + ", ".join(c_table(i) for i in table) + "}"


def tables_to_c(tables: dict[str, np.ndarray]) -> str:
    """Convert tables of values to a string of code."""
    return "\n".join(
        f"static const double {variable}["
        + "][".join(f"{i}" for i in table.shape)
        + "] = "
        + c_table(table)
        + ";"
        for variable, table in tables.items()
    )


def generate(
    form: RepresentedByGraph,
    language: str = "C",
) -> tuple[str, dict[RepresentedByGraph, str]]:
    """Generate code.

    Args:
        form: The form or other object to be assembled
        language: The programming language to use

    Returns:
        Code
    """
    if language != "C":
        raise NotImplementedError("Only generation of C is supported for now")

    graph = form.graph

    assert is_dag(graph)

    # TODO: get this from somewhere
    rules = {
        dx: QuadratureRule([[1 / 6, 1 / 6], [2 / 3, 1 / 6], [1 / 6, 2 / 3]], [1 / 6, 1 / 6, 1 / 6])
    }

    tables = {}

    print("Graph:")
    print_graph(graph)
    print()

    print("Applying integrals_to_quadrature")
    graph = integrals_to_quadrature(graph, rules)
    print()

    print("Graph:")
    print_graph(graph)
    print()

    print("Applying apply_push_forwards")
    graph = apply_push_forwards(graph)
    print()

    print("Graph:")
    print_graph(graph)
    print()

    print("Applying tabulate_finite_elements")
    q_tables, graph = tabulate_quadrature(graph)
    tables = {**tables, **q_tables}
    print()

    print("Tables:")
    print(tables)
    print("Graph:")
    print_graph(graph)
    print()

    print("Applying expand_jacobians")
    graph = expand_jacobians(graph)
    print()

    print("Graph:")
    print_graph(graph)
    print()

    print("Applying tabulate_finite_elements")
    fe_tables, graph = tabulate_finite_elements(graph)
    tables = {**tables, **fe_tables}
    print()

    print("Tables:")
    print(tables)
    print("Graph:")
    print_graph(graph)
    print()

    print("Applying convert_complex_to_real")
    graph = convert_complex_to_real(graph)
    print()

    print("Graph:")
    print_graph(graph)
    print()

    code = (
        "void tabulate_tensor_f64(\n"
        "    double* restrict A, const double* restrict w, const double* restrict c,\n"
        "    const double* restrict coordinate_dofs,\n"
        "    const int* restrict entity_local_index,\n"
        "    const uint8_t* restrict quadrature_permutation, void* custom_data\n"
        ") {\n"
    )

    code += indented(tables_to_c(tables), 2)
    code += "\n\n"
    code += indented(graph.root.generate_c(), 2)
    code += "\n}\n"

    """
    code += "  static const double weights[3] = {"
    code += ", ".join(f"{w}" for w in rules[dx].weights)
    code += "};\n"
    code += (
        f"  static const double basis_values[1][1][{table.shape[0]}][{table.shape[1]}] = {{{{{{\n"
    )
    code += ",\n".join("    {" + ", ".join(f"{i}" for i in row) + "}" for row in table)
    code += "}}};"
    code += (
        "  static const double FE1_C0_D10_Q48e[1][1][1][3] = {{{{-1.0, 1.0, 0.0}}}};\n"
        "  static const double FE1_C1_D01_Q48e[1][1][1][3] = {{{{-1.0, 0.0, 1.0}}}};\n"
        "  double J0_c0 = 0.0;\n"
        "  double J0_c3 = 0.0;\n"
        "  double J0_c1 = 0.0;\n"
        "  double J0_c2 = 0.0;\n"
        "  for (int ic = 0; ic < 3; ++ic)\n"
        "  {\n"
        "    J0_c0 += coordinate_dofs[ic * 3] * FE1_C0_D10_Q48e[0][0][0][ic];\n"
        "    J0_c3 += coordinate_dofs[ic * 3 + 1] * FE1_C1_D01_Q48e[0][0][0][ic];\n"
        "    J0_c1 += coordinate_dofs[ic * 3] * FE1_C1_D01_Q48e[0][0][0][ic];\n"
        "    J0_c2 += coordinate_dofs[ic * 3 + 1] * FE1_C0_D10_Q48e[0][0][0][ic];\n"
        "  }\n"
        "  double jdet = J0_c0 * J0_c3 - J0_c1 * J0_c2;\n"
        "  double abs_jdet = fabs(jdet);\n"
        "  for (int iq = 0; iq < 3; ++iq)\n"
        "  {\n"
        "    double fw0 = abs_jdet * weights[iq];\n"
        "    double temp_0[3] = {0};\n"
        "    for (int j = 0; j < 3; ++j)\n"
        "      temp_0[j] = fw0 * basis_values[0][0][iq][j];\n"
        "    for (int j = 0; j < 3; ++j)\n"
        "      for (int i = 0; i < 3; ++i)\n"
        "        A[3 * i + j] += basis_values[0][0][iq][i] * temp_0[j];\n"
        "  }\n"
        "}\n"
    )
    """
    signatures = {
        form: (
            "void tabulate_tensor_f64(double* restrict, const double* restrict, "
            "const double* restrict, const double* restrict, const int* restrict, "
            "const uint8_t* restrict, void*);"
        ),
    }

    return code, signatures
