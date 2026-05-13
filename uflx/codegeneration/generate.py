"""Code generation."""
from typing import Self, Protocol
import numpy as np
from uflx.integrals import AbstractIntegral, dx
from uflx.graphs import RepresentedByGraph, is_dag, Graph, GraphNode, generate_graph
from uflx.graphs.algorithms import replace
from uflx.functions import Argument
from uflx.function_spaces import AbstractReferenceMappedFunctionSpace
from uflx.operators import Conj, Mult
import networkx as nx


class ConvertToCode(Protocol):
    def as_code(self, language: str, bracketed: bool = False) -> str:
        """Generate code for this object."""


def reconstruct(object, args, replacements):
    if all(a not in replacements for a in args):
        return object
    return object.__class__(*(replacements.get(a, a) for a in args))


def print_node(graph: Graph, node: GraphNode, indentation: int = 0):
    print(" " * (2 * indentation) + f"{node!r}")
    for next in graph.successors(node):
        print_node(graph, next, indentation + 1)


def print_graph(graph: Graph):
    print_node(graph, graph.root)


class QuadratureRule:
    def __init__(self, points, weights):
        self.points = points
        self.weights = weights

    @property
    def npoints(self):
        return len(self.weights)


class QuadraturePoint:
    def __init__(self, rule, index):
        self.rule = rule
        self.index = index

    def __repr__(self):
        return f"QuadraturePoint({self.index})"


class QuadratureWeight:
    def __init__(self, rule, index):
        self.rule = rule
        self.index = index

    def __repr__(self):
        return f"QuadratureWeight({self.index})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self


class AbsJacobianDeterminant:
    def __init__(self, domain):
        self.domain = domain

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return self

    def as_code(self, language: str, bracketed: bool = False) -> str:
        """Generate code for this object."""
        match language:
            case "C":
                return "TODO"
        raise NotImplementedError()


class QuadratureLoop:
    def __init__(self, body, rule, variable):
        self.body = body
        self.rule = rule
        self.variable = variable

    def __repr__(self):
        return f"QuadratureLoop({self.variable})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return reconstruct(self, (self.body, self.rule, self.variable), replacements)

    def as_code(self, language: str, bracketed: bool = False) -> str:
        """Generate code for this object."""
        match language:
            case "C":
                return (
                    f"for (int {self.variable}=0; {self.variable}!={self.rule.npoints}; {self.variable}++)\n"
                    "{\n" + indented(self.body.as_code(language), 2) + "\n}"
                )
        raise NotImplementedError()


def indented(code, spaces):
    return "\n".join(" " * spaces + line for line in code.split("\n"))


class Loop:
    def __init__(self, variable: str, start: int, end: int, body):
        self.variable = variable
        self.start = start
        self.end = end
        self.body = body

    def __repr__(self):
        return f"Loop({self.variable}, {self.start}, {self.end})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.body}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return reconstruct(self, (self.variable, self.start, self.end, self.body), replacements)

    def as_code(self, language: str, bracketed: bool = False) -> list[str]:
        """Generate code for this object."""
        match language:
            case "C":
                return (
                    f"for (int {self.variable}={self.start}; {self.variable}!={self.end};{self.variable}++)\n"
                    "{\n" + indented(self.body.as_code(language), 2) + "\n}"
                )
        raise NotImplementedError()


class EvaluatedBasisFunction:
    def __init__(self, element, basis_index, point):
        self.element = element
        self.basis_index = basis_index
        self.point = point

    def __repr__(self):
        return f"EvaluatedBasisFunction({self.element!r}, {self.basis_index}, {self.point!r})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return reconstruct(self, (self.element, self.basis_index, self.point), replacements)


class PushedForward:
    def __init__(self, map, function):
        self.map = map
        self.function = function

    def __repr__(self):
        return f"PushedForward({self.map})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.function}

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""
        return reconstruct(self, (self.map, self.function), replacements)


class VariableNamer:
    def __init__(self):
        self.i = -1
        self.n = -1
        self.fe_i = -1
        self.qr_i = -1

    def variable(self):
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
        self.fe_i += 1
        return f"FE{self.fe_i}"

    def quadrature_table(self):
        self.qr_i += 1
        return f"QW{self.qr_i}"


global_variable_namer = VariableNamer()
local_tensor = "A"


def flatten_component(indices, shape, bracketed=False):
    assert len(indices) == len(shape)
    if len(indices) == 1:
        return indices[0]

    component = flatten_component(indices[:-1], shape[:-1], True) + f" * {shape[-1]} + {indices[-1]}"
    if bracketed:
        return f"({component})"
    else:
        return f"{component}"


class AddToLocalTensor:
    def __init__(self, component, shape, body):
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
        return f"AddToLocalTensor({self.component})"

    def as_code(self, language: str, bracketed: bool = False) -> str:
        """Generate code for this object."""
        match language:
            case "C":
                return f"{local_tensor}[" + flatten_component(self.component, self.shape) + "] += " + self.body.as_code(language) + ";"
        raise NotImplementedError()


class ArrayEntry:
    def __init__(self, array, index):
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
        return f"ArrayEntry({self.array}, {self.index})"

    def as_code(self, language: str, bracketed: bool = False) -> str:
        """Generate code for this object."""
        match language:
            case "C":
                return f"{self.array}[" + "][".join(self.index) + "]"
        raise NotImplementedError()


def apply_push_forwards(
    graph: Graph,
) -> Graph:
    return replace(graph, {node: node.map.push_forward_symbolic(node.function) for node in graph if isinstance(node, PushedForward)})


def convert_complex_to_real(
    graph: Graph,
) -> Graph:
    return replace(graph, {node: node.init_arg for node in graph if isinstance(node, Conj)})


def integrals_to_quadrature(
    graph: Graph,
    rules: dict[RepresentedByGraph, QuadratureRule],
    variable_namer=global_variable_namer,
) -> Graph:
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
                    raise NotImplementedError("Code generation only implemented for reference mapped spaces")
                if len(i_space.elements) != 1:
                    raise NotImplementedError("Code generation currently only implemented for spaces with exactly one element")
                tensor_shape_components[i.component] = i_space.elements[0].dim
            tensor_shape = tuple(
                tensor_shape_components.get(i, 1)
                for i in range(min(tensor_shape_components.keys()), max(tensor_shape_components.keys()) + 1)
            )
            variables = tuple(0 if component == 1 else variable_namer.variable() for component in tensor_shape)

            qvariable = variable_namer.variable()

            for a in arguments:
                to_replace[a] = PushedForward(a.function_space.elements[0].map_type, EvaluatedBasisFunction(a.function_space.elements[0], qvariable, QuadraturePoint(rules[node.measure], variables[a.component])))

            domain = arguments[0].function_space.domain
            for a in arguments:
                assert a.function_space.domain == domain

            integrand = Mult(Mult(QuadratureWeight(rules[node.measure], variables[a.component]), AbsJacobianDeterminant(domain)), node.integrand)

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
    table_map = {}
    tables = {}
    to_replace = {}
    for node in graph:
        if isinstance(node, EvaluatedBasisFunction):
            if not isinstance(node.point, QuadraturePoint):
                raise NotImplementedError()

            id = (node.element, node.point.rule)
            if id not in table_map:
                name = variable_namer.finite_element_table()
                table_map[id] = name
                tables[name] = node.element.tabulate(node.point.rule.points)
            to_replace[node] = ArrayEntry(table_map[id], (node.point.index, node.basis_index))

    return tables, replace(graph, to_replace)


def tabulate_quadrature(
    graph,
    variable_namer=global_variable_namer,
):
    table_map = {}
    tables = {}
    to_replace = {}
    for node in graph:
        if isinstance(node, QuadratureWeight):
            if node.rule not in table_map:
                name = variable_namer.quadrature_table()
                table_map[node.rule] = name
                tables[name] = np.array(node.rule.weights)
            to_replace[node] = ArrayEntry(table_map[node.rule], (node.index, ))

    return tables, replace(graph, to_replace)


def c_table(table):
    if len(table.shape) == 1:
        return "{" + ", ".join(f"{i}" for i in table) + "}"
    return "{" + ", ".join(c_table(i) for i in table) + "}"

def tables_to_code(tables, language: str):
    match language:
        case "C":
            return "\n".join(
                f"static const double {variable}[" + "][".join(f"{i}" for i in table.shape) + "] = " + c_table(table) + ";"
                for variable, table in tables.items()
            )
    raise NotImplementedError()


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
    qtables, graph = tabulate_quadrature(graph)
    print()

    print("Tables:")
    print(qtables)
    print("Graph:")
    print_graph(graph)
    print()

    print("Applying tabulate_finite_elements")
    tables, graph = tabulate_finite_elements(graph)
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

    print(tables_to_code({**tables, **qtables}, language))

    print(graph.root.as_code(language))


    # TODO: get these from Basix
    table = np.array([[1 - x - y, x, y] for x, y in rules[dx].points])
    table_deriv_0 = np.array([-1.0, 1.0, 0.0])
    table_deriv_1 = np.array([-1.0, 0.0, 1.0])

    code = (
        "void tabulate_tensor_f64(\n"
        "    double* restrict A, const double* restrict w, const double* restrict c,\n"
        "    const double* restrict coordinate_dofs,\n"
        "    const int* restrict entity_local_index,\n"
        "    const uint8_t* restrict quadrature_permutation, void* custom_data\n"
        ") {\n"
    )
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

    signatures = {
        form: (
            "void tabulate_tensor_f64(double* restrict, const double* restrict, "
            "const double* restrict, const double* restrict, const int* restrict, "
            "const uint8_t* restrict, void*);"
        ),
    }

    return code, signatures
