"""Code generation."""

import numpy as np
from uflx.integrals import AbstractIntegral, dx
from uflx.graphs import RepresentedByGraph, is_dag, Graph, GraphNode
from uflx.graphs.algorithms import replace
from uflx.functions import Argument
from uflx.function_spaces import AbstractReferenceMappedFunctionSpace
import networkx as nx


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


class QuadraturePoint:
    def __init__(self, rule, index):
        self.rule = rule
        self.index = index

    def __repr__(self):
        return f"QuadraturePoint({self.index})"

class QuadratureLoop:
    def __init__(self, integrand, rule, variable):
        self.integrand = integrand
        self.rule = rule
        self.variable = variable

    def __repr__(self):
        return f"QuadratureLoop({self.variable})"


class Loop:
    def __init__(self, variable: str, start: int, end: int, body):
        self.variable = variable
        self.start = start
        self.end = end

    def __repr__(self):
        return f"Loop({self.variable}, {self.start}, {self.end})"


class EvaluatedBasisFunction:
    def __init__(self, element, basis_index, point):
        self.element = element
        self.basis_index = basis_index
        self.point = point

    def __repr__(self):
        return f"EvaluatedBasisFunction({self.element!r}, {self.basis_index}, {self.point!r})"


class PushedForward:
    def __init__(self, function):
        self.function = function


    def __repr__(self):
        return f"PushedForward({self.function!r})"


def integrals_to_quadrature(
    graph: Graph,
    rules: dict[RepresentedByGraph, QuadratureRule],
) -> Graph:

    print("Graph of integral:")
    print_graph(graph)
    print()

    new_graph = Graph()

    updated_nodes = {}
    to_update = {}
    for node in reversed(list(nx.topological_sort(graph))):
        if isinstance(node, AbstractIntegral):
            tensor_shape_components = {}

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
            assert len(tensor_shape) <= 6
            variables = ["i", "j", "k", "l", "m", "n"][:len(tensor_shape)]

            qvariable = "qi"

            for a in arguments:
                to_update[a] = PushedForward(EvaluatedBasisFunction(a.function_space.elements[0], qvariable, QuadraturePoint(rules[node.measure], variables[a.component])))

            from IPython import embed; embed()


            qloop = QuadratureLoop(node.integrand, rules[node.measure], qvariable)
            if graph.root == node:
                new_graph.set_root(qloop)
            new_graph.add_node(qloop)
            for pre in graph.predecessors(node):
                new_graph.add_edge(updated_nodes.get(pre, pre), qloop)
            for post in graph.successors(node):
                if post != node.measure:
                    new_graph.add_edge(qloop, updated_nodes.get(post, post))

            next = qloop
            for variable, count in zip(variables[::-1], tensor_shape[::-1]):
                loop = Loop(variable, 0, count, next)
                new_graph.add_node(loop)
                new_graph.add_edge(loop, next)
                if new_graph.has_root and new_graph.root == next:
                    new_graph.set_root(loop)
                next = loop

            updated_nodes[node] = next
        else:
            new_graph.add_node(node)
            for post in graph.successors(node):
                new_graph.add_edge(node, updated_nodes.get(post, post))
            if graph.root == node:
                new_graph.set_root(node)

    # TODO: add child nodes too inside replace
    new_graph = replace(new_graph, to_update)

    print("Graph of computation:")
    print_graph(new_graph)
    print()

    # from IPython import embed; embed()

    return new_graph


def generate(
    form: RepresentedByGraph, language: str = "C"
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

    graph = integrals_to_quadrature(graph, rules)

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
