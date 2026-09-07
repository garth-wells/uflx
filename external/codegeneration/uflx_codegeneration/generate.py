"""Code generation."""

import numpy as np
import numpy.typing as npt
import quadraturerules
from uflx.algorithms import pull_back_to_reference, replace
from uflx.geometry import (
    expand_geometry,
)
from uflx.graphs import (
    GraphNode,
    as_graph,
)
from uflx.integrals import AbstractMeasure, dx
from uflx.maps import apply_push_forwards
from uflx.points import Point

from uflx_codegeneration import symbols
from uflx_codegeneration.algorithms import (
    expand_inner_products,
    insert_geometry_functions,
    tabulate_finite_elements,
)
from uflx_codegeneration.c import GenerateC, tables_to_c
from uflx_codegeneration.nodes import ArrayEntry
from uflx_codegeneration.quadrature import (
    QuadraturePoint,
    QuadratureRule,
    QuadratureWeight,
    integrals_to_quadrature,
    quadrature_rule,
)
from uflx_codegeneration.utils import indented


def tabulate_quadrature(
    expression: GraphNode,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, npt.NDArray(np.floating)], GraphNode]:
    """Generate tables of values for quadrature rules."""
    table_map = {}
    tables = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    for node in as_graph(expression):
        if isinstance(node, QuadratureWeight):
            id = (node.rule, "weights")
            if id not in table_map:
                name = variable_namer.quadrature_table()
                table_map[id] = name
                tables[name] = node.rule.weights
            to_replace[node] = ArrayEntry(table_map[id], (node.index,))
        if isinstance(node, QuadraturePoint):
            id = (node.rule, "points")
            if id not in table_map:
                name = variable_namer.quadrature_table()
                table_map[id] = name
                tables[name] = node.rule.points
            to_replace[node] = Point(
                [
                    ArrayEntry(
                        table_map[id],
                        (
                            node.dim * node.index + i
                            if isinstance(node.index, int)
                            else f"{node.dim} * {node.index} + {i}",
                        ),
                    )
                    for i in range(node.dim)
                ]
            )

    return tables, replace(expression, to_replace)


def generate(
    form: GraphNode,
    language: str = "C",
) -> tuple[str, dict[GraphNode, str]]:
    """Generate code.

    Args:
        form: The form or other object to be assembled
        language: The programming language to use

    Returns:
        Code
    """
    if language != "C":
        raise NotImplementedError("Only generation of C is supported for now")

    # TODO: get this from somewhere
    rules: dict[AbstractMeasure, QuadratureRule] = {}
    # For now, use a degree 10 rule:
    points, weights = quadraturerules.single_integral_quadrature(
        quadraturerules.QuadratureRule.XiaoGimbutas,
        quadraturerules.Domain.Triangle,
        10,
    )
    rules[dx] = quadrature_rule([p[1:] for p in points], 0.5 * weights)

    # Apply algorithms from UFLx
    form = pull_back_to_reference(form)
    form = apply_push_forwards(form)

    # Apply codegeneration algorithms
    form = integrals_to_quadrature(form, rules)
    geometry_functions, form = insert_geometry_functions(form)
    form = expand_geometry(form)
    form = expand_inner_products(form)

    # Tabulate quadrature rules and finite element functions
    q_tables, form = tabulate_quadrature(form)
    fe_tables, form = tabulate_finite_elements(form)
    tables = {**q_tables, **fe_tables}

    code = ""
    for fname, (dtype, inputs, function) in geometry_functions.items():
        code += f"{dtype} {fname}("
        code += ", ".join(f"{i._dtype} {i._variable}" for i in inputs)
        code += ") {\n"
        ftables, function = tabulate_finite_elements(function)
        code += indented(tables_to_c(ftables), 2)
        code += "\n\n"
        assert isinstance(function, GenerateC)
        code += f"  return {function.generate_c()};\n"
        code += "}\n\n"
    code += (
        "void tabulate_tensor_f64(\n"
        f"    double* restrict {symbols.local_tensor},\n"
        f"    const double* restrict {symbols.coefficients},\n"
        f"    const double* restrict {symbols.constants},\n"
        f"    const double* restrict {symbols.coordinate_dofs},\n"
        f"    const int* restrict {symbols.entity_local_index},\n"
        f"    const uint8_t* restrict {symbols.quadrature_permutation},\n"
        f"    void* {symbols.custom_data}\n"
        ") {\n"
    )

    code += indented(tables_to_c(tables), 2)
    code += "\n\n"
    assert isinstance(form, GenerateC)
    code += indented(form.generate_c(), 2)
    code += "\n}\n"

    signatures = {
        form: (
            "void tabulate_tensor_f64(double* restrict, const double* restrict, "
            "const double* restrict, const double* restrict, const int* restrict, "
            "const uint8_t* restrict, void*);"
        ),
    }

    return code, signatures
