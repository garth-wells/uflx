"""Code generation."""

import quadraturerules
from uflx.basis_functions import EvaluatedPhysicalBasisFunction, EvaluatedReferenceBasisFunction
from uflx.complex import take_real_part
from uflx.domains import AbstractCoordinateElement, AbstractDomain
from uflx.function_spaces import AbstractReferenceMappedFunctionSpace
from uflx.functions import AbstractPhysicalFunction, Argument, ReferenceArgument
from uflx.geometry import (
    JacobianDeterminant,
    JacobianInverseTranspose,
    JacobianInverse,
    JacobianTranspose,
    Jacobian,
    PhysicalToReference,
    ReferenceToPhysical,
    SingleSpatialCoordinate,
    expand_geometry,
)
from uflx.graphs import (
    Graph,
    GraphNode,
    RepresentedByGraph,
    generate_graph,
)
from uflx.finite_elements import AbstractReferenceMappedFiniteElement
from uflx.graphs.algorithms import reconstruct_node, replace, pull_back_to_reference
from uflx.integrals import AbstractIntegral, AbstractMeasure, Measure, dx
from uflx.maps import PushedForward, apply_push_forwards
from uflx.operators import Grad, ReferenceGrad
from uflx.points import Point, PointComponent

from uflx_codegeneration import symbols
from uflx_codegeneration.algorithms import (
    expand_inner_products,
    insert_geometry_functions,
    tabulate_finite_elements,
)
from uflx_codegeneration.c import GenerateC, tables_to_c
from uflx_codegeneration.nodes import AddToLocalTensor, ArrayEntry, Loop
from uflx_codegeneration.quadrature import (
    QuadratureLoop,
    QuadraturePoint,
    QuadratureRule,
    QuadratureWeight,
    quadrature_rule,
)
from uflx_codegeneration.utils import indented


def extract_domain(graph: Graph, node: GraphNode) -> AbstractDomain:
    """Extract the domain associated with a node."""
    domain: AbstractDomain | None = None
    for i in graph.descendants(node):
        if isinstance(i, AbstractPhysicalFunction):
            if domain is None:
                domain = i.function_space.domain
            else:
                assert domain == i.function_space.domain
        if hasattr(i, "domain"):
            if domain is None:
                domain = i.domain
            else:
                assert domain == i.domain
    assert domain is not None
    return domain


def integrals_to_quadrature(
    graph: Graph,
    rules: dict[AbstractMeasure, QuadratureRule],
    variable_namer=symbols.global_variable_namer,
) -> Graph:
    """Replace integrals with quadrature."""
    updated_nodes: dict[GraphNode, GraphNode] = {}
    to_replace: dict[GraphNode, GraphNode] = {}

    for node in graph.ordered_nodes():
        if isinstance(node, AbstractIntegral):
            rule = rules[node.measure]
            qvariable = variable_namer.variable()
            qpoint = QuadraturePoint(rule, qvariable)

            tensor_shape_components = {}

            if not isinstance(node.measure, Measure):
                raise NotImplementedError()
            if node.measure._codim != 0 or node.measure._boundary_only:
                raise NotImplementedError("Only codim 0 integrals supported for now")

            arguments = []
            for i in graph.descendants(node):
                if isinstance(i, (Argument, ReferenceArgument)) and i.integral_label == node.label:
                    arguments.append(i)
                if isinstance(i, SingleSpatialCoordinate):
                    domain = extract_domain(graph, node)
                    if not isinstance(domain, AbstractCoordinateElement):
                        raise NotImplementedError(
                            "Code generation only implemented for reference mapped domain"
                        )
                    to_replace[i] = PointComponent(ReferenceToPhysical(qpoint, domain), i._component)
                if isinstance(i, Jacobian) and i.point is None:
                    to_replace[i] = Jacobian(i.domain, qpoint)
                if isinstance(i, JacobianInverse) and i.point is None:
                    to_replace[i] = JacobianInverse(i.domain, qpoint)
                if isinstance(i, JacobianTranspose) and i.point is None:
                    to_replace[i] = JacobianTranspose(i.domain, qpoint)
                if isinstance(i, JacobianInverseTranspose) and i.point is None:
                    to_replace[i] = JacobianInverseTranspose(i.domain, qpoint)
                if isinstance(i, JacobianDeterminant) and i.point is None:
                    to_replace[i] = JacobianDeterminant(i.domain, qpoint)
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
                tensor_shape_components[i.component_index] = i_space.elements[0].dim
            tensor_shape = tuple(
                tensor_shape_components.get(i, 1)
                for i in range(
                    min(tensor_shape_components.keys()), max(tensor_shape_components.keys()) + 1
                )
            )
            variables = tuple(
                "0" if component == 1 else variable_namer.variable() for component in tensor_shape
            )

            for a in arguments:
                assert isinstance(a.function_space, AbstractReferenceMappedFunctionSpace)
                assert isinstance(a.function_space.domain, AbstractCoordinateElement)
                if isinstance(a, Argument):
                    to_replace[a] = EvaluatedPhysicalBasisFunction(
                        a.function_space,
                        a.function_space.elements[0],
                        variables[a.component_index],
                        ReferenceToPhysical(qpoint, a.function_space.domain),
                    )
                elif isinstance(a, ReferenceArgument):
                    to_replace[a] = EvaluatedReferenceBasisFunction(
                        a.function_space.elements[0],
                        variables[a.component_index],
                        qpoint,
                    )

            domain = arguments[0].function_space.domain
            for a in arguments:
                assert a.function_space.domain == domain

            assert isinstance(domain, AbstractCoordinateElement)
            integrand = (
                QuadratureWeight(rules[node.measure], qvariable)
                # * abs(JacobianDeterminant(domain, QuadraturePoint(rules[node.measure], qvariable)))
                * node.integrand
            )

            body = AddToLocalTensor(variables, tensor_shape, integrand)

            qloop = QuadratureLoop(body, rules[node.measure], qvariable)

            next: GraphNode = qloop
            for variable, count in zip(variables[::-1], tensor_shape[::-1]):
                if variable == "0":
                    continue
                assert isinstance(count, int)
                loop = Loop(variable, 0, count, next)
                next = loop

            updated_nodes[node] = next

    new_graph = generate_graph(updated_nodes.get(graph.root, graph.root))
    new_graph = replace(new_graph, to_replace)
    return new_graph


def tabulate_quadrature(
    graph,
    variable_namer=symbols.global_variable_namer,
):
    """Generate tables of values for quadrature rules."""
    table_map = {}
    tables = {}
    to_replace: dict[GraphNode, GraphNode] = {}
    for node in graph:
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

    return tables, replace(graph, to_replace)


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

    assert graph.is_dag()

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
    graph = pull_back_to_reference(graph)
    graph = apply_push_forwards(graph)

    # Apply codegeneration algorithms
    graph = integrals_to_quadrature(graph, rules)
    geometry_functions, graph = insert_geometry_functions(graph)
    graph = expand_geometry(graph)
    graph = expand_inner_products(graph)
    graph = take_real_part(graph)

    q_tables, graph = tabulate_quadrature(graph)
    fe_tables, graph = tabulate_finite_elements(graph)
    tables = {**q_tables, **fe_tables}

    code = ""
    for fname, (dtype, inputs, fgraph) in geometry_functions.items():
        code += f"{dtype} {fname}("
        code += ", ".join(f"{i._dtype} {i._variable}" for i in inputs)
        code += ") {\n"
        ftables, fgraph = tabulate_finite_elements(fgraph)
        code += indented(tables_to_c(ftables), 2)
        code += "\n\n"
        assert isinstance(fgraph.root, GenerateC)
        code += f"  return {fgraph.root.generate_c()};\n"
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
    assert isinstance(graph.root, GenerateC)
    code += indented(graph.root.generate_c(), 2)
    code += "\n}\n"

    signatures = {
        form: (
            "void tabulate_tensor_f64(double* restrict, const double* restrict, "
            "const double* restrict, const double* restrict, const int* restrict, "
            "const uint8_t* restrict, void*);"
        ),
    }

    return code, signatures
