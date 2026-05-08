"""Code generation."""

import numpy as np
from uflx.integrals import AbstractIntegral, dx
from uflx.graphs import RepresentedByGraph, is_dag, Graph
from uflx.functions import FiniteElementFunction


class QuadratureRule:
    def __init__(self, points, weights):
        self.points = points
        self.weights = weights


class QuadratureLoop:
    def __init__(self, integrand, quadrature_rule):
        self.integrand = integrand
        self.quadrature_rule = quadrature_rule


def integrals_to_quadrature(
    graph: Graph,
    rules: dict[RepresentedByGraph, QuadratureRule],
) -> Graph:

    new_graph = Graph(graph)

    updated_nodes = {}
    for node in graph:
        if isinstance(node, AbstractIntegral):
            new_graph.remove_node(node)
            if len(list(new_graph.predecessors(node.measure))) == 0:
                new_graph.remove_node(node.measure)
            loop = QuadratureLoop(node.integrand, rules[node.measure])
            new_graph.add_node(loop)
            for pre in graph.predecessors(node):
                new_graph.add_edge(updated_nodes.get(pre, pre), loop)
            for post in graph.successors(node):
                if post != node.measure:
                    new_graph.add_edge(loop, updated_nodes.get(post, post))

            tensor_shape = {}

            # for i in networkx.descendents(graph, node):
            #    if isinstance(i, FiniteElementFunction):
            #        pass

            updated_nodes[node] = loop

            print(node)
            print(rules[node.measure])
            print(list(graph.nodes))
            print(list(node.graph.nodes))

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
