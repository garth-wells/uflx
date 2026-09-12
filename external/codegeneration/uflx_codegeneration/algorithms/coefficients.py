"""Coefficient algorithms."""

from typing import Any

from uflx.algorithms import replace
from uflx.basis_functions import EvaluatedReferenceBasisFunction
from uflx.graphs import GraphNode, as_graph

from uflx_codegeneration import symbols
from uflx_codegeneration.coefficients import EvaluatedReferenceCoefficientBasisFunction
from uflx_codegeneration.nodes import (
    AccumulateToVariable,
    ArrayEntry,
    Block,
    Declare,
    FunctionCall,
    Loop,
    Return,
    Variable,
)


def insert_coefficient_functions(
    expression: GraphNode,
    variable_namer: symbols.VariableNamer = symbols.global_variable_namer,
) -> tuple[dict[str, tuple[str, list[Variable], GraphNode]], GraphNode]:
    """Replace evaluated-coefficient basis functions with calls to summation functions.

    Each distinct EvaluatedReferenceCoefficientBasisFunction found in the
    expression (one per distinct combination of coefficient, derivative and
    component -- see that class) is replaced with a call to a small,
    self-contained C function that sums the coefficient's degrees of freedom
    against the appropriate finite element table, eg:

        static double coeff0(const double* restrict w, int q0) {
          static const double FE0[...][...][...] = {...};
          double acc0 = 0.0;
          for (int k0=0; k0!=3; ++k0) { acc0 += (FE0[0][q0][k0] * w[k0]); }
          return acc0;
        }

    This mirrors ``insert_geometry_functions``: the returned dict maps
    function names to (return type, input variables, function body), and it
    is the caller's job (see ``generate.generate``) to tabulate any finite
    elements referenced in each function body (via ``tabulate_finite_elements``)
    and emit the function definitions before they are used.

    This must be called after ``expand_geometry`` and ``expand_inner_products``
    have fully resolved any derivatives/components of the coefficient (ie
    called ``.diff()``/``.component()`` on it), since only then is it known
    which distinct (coefficient, derivative, component) combinations actually
    need a summation loop.
    """
    functions: dict[str, tuple[str, list[Variable], GraphNode]] = {}
    to_replace: dict[GraphNode, GraphNode] = {}

    nodes = [
        node
        for node in as_graph(expression)
        if isinstance(node, EvaluatedReferenceCoefficientBasisFunction)
    ]

    # Every distinct coefficient (identified by its count) is given a
    # contiguous block of `ndofs` slots in the coefficients array, ordered by
    # count. (There is currently no other convention establishing the layout
    # of the coefficients array, since this is the first code generation
    # support for Coefficients.)
    ndofs_by_count: dict[int, int] = {}
    for node in nodes:
        ndofs_by_count[node.count] = node.element.dim
    offset_by_count: dict[int, int] = {}
    offset = 0
    for count in sorted(ndofs_by_count):
        offset_by_count[count] = offset
        offset += ndofs_by_count[count]

    for node in nodes:
        fname = variable_namer.coefficient_function_name()
        dof_variable = node.basis_index
        assert isinstance(dof_variable, str)
        accumulator = variable_namer.variable()

        table_lookup = EvaluatedReferenceBasisFunction(
            node.element,
            dof_variable,
            node.point,
            node.derivative,
            node.component_index,
        )
        dof_offset = offset_by_count[node.count]
        w_index = dof_variable if dof_offset == 0 else f"{dof_offset} + {dof_variable}"

        w = Variable("const double* restrict", symbols.coefficients)
        inputs: list[Variable] = [w]
        f_args: list[Any] = [symbols.coefficients]
        if isinstance(node.point_index, str):
            inputs.append(Variable("int", node.point_index))
            f_args.append(node.point_index)

        body = Block(
            (
                Declare("double", accumulator, 0.0),
                Loop(
                    dof_variable,
                    0,
                    node.element.dim,
                    AccumulateToVariable(
                        accumulator,
                        table_lookup * ArrayEntry(symbols.coefficients, (w_index,)),
                    ),
                ),
                Return(Variable("double", accumulator)),
            )
        )
        functions[fname] = ("double", inputs, body)
        to_replace[node] = FunctionCall(fname, *f_args)

    return functions, replace(expression, to_replace)
