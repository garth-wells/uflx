"""Experimental extraction of cellwise geometry from lowered UFLx graphs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import basix
from uflx.expressions import (
    Abs,
    AbstractExpression,
    AbstractScalar,
    MatVec,
    Mult,
    expression_sum,
)
from uflx.geometry import JacobianDeterminant, JacobianInverseTranspose
from uflx.graphs import Graph, GraphNode
from uflx.graphs.algorithms import replace
from uflx.operators import Inner, ReferenceGrad


@dataclass(frozen=True)
class GeometryKernelSpec:
    """Description of a geometry-only kernel required by a lowered form."""

    output_size: int = 6
    scope: str = "cell"


class GeometryTensorComponent(AbstractScalar):
    """One component of a symmetric 3-by-3 cellwise geometry tensor.

    The tensor is ``G = abs(det(J)) * inv(J) * inv(J).T`` and is stored in
    packed upper-triangular order ``(G00, G01, G02, G11, G12, G22)``.
    """

    def __init__(self, row: int, column: int):
        """Initialise a packed symmetric tensor component."""
        if not 0 <= row < 3 or not 0 <= column < 3:
            raise IndexError("geometry tensor component index out of range")
        self.row = min(row, column)
        self.column = max(row, column)

    @property
    def packed_index(self) -> int:
        """Return the component's upper-triangular packed index."""
        return {(0, 0): 0, (0, 1): 1, (0, 2): 2, (1, 1): 3, (1, 2): 4, (2, 2): 5}[
            (self.row, self.column)
        ]

    @property
    def successors(self) -> set[GraphNode]:
        """Return this leaf node's empty successor set."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """Return arguments needed to reconstruct this node."""
        return self.row, self.column

    def __repr__(self) -> str:
        """Return a compact representation."""
        return f"GeometryTensorComponent({self.row}, {self.column})"


def geometry_kernel_name(kernel_name: str) -> str:
    """Return the exported geometry-kernel symbol paired with ``kernel_name``."""
    return f"{kernel_name}_geometry"


def _poisson_metric_contraction(node: GraphNode):
    """Return a reference-gradient contraction when ``node`` is affine Poisson geometry."""
    if not isinstance(node, Mult):
        return None

    if isinstance(node.first, Abs) and isinstance(node.first.argument, JacobianDeterminant):
        determinant, contraction = node.first.argument, node.second
    elif isinstance(node.second, Abs) and isinstance(node.second.argument, JacobianDeterminant):
        determinant, contraction = node.second.argument, node.first
    else:
        return None

    if not isinstance(contraction, Inner):
        return None
    if not isinstance(contraction.first, MatVec) or not isinstance(contraction.second, MatVec):
        return None

    left, right = contraction.first, contraction.second
    if not isinstance(left.first, JacobianInverseTranspose) or not isinstance(
        right.first, JacobianInverseTranspose
    ):
        return None
    if not isinstance(left.second, ReferenceGrad) or not isinstance(right.second, ReferenceGrad):
        return None
    if left.first.domain is not determinant.domain or right.first.domain is not determinant.domain:
        return None
    if left.first.value_shape != (3, 3) or right.first.value_shape != (3, 3):
        return None

    domain = determinant.domain
    if len(domain.elements) != 1 or domain.elements[0].lagrange_superdegree != 1:
        return None

    grad_left = left.second.expand_geometry()
    grad_right = right.second.expand_geometry()
    assert isinstance(grad_left, AbstractExpression)
    assert isinstance(grad_right, AbstractExpression)
    # Preserve the matrix-vector factorisation instead of flattening this
    # immediately into nine scalar terms. The code generator can then
    # precompute the three components of G @ grad_right once per trial dof,
    # rather than materialising one scratch vector per tensor entry.
    return expression_sum(
        grad_left.component(i)
        * expression_sum(GeometryTensorComponent(i, j) * grad_right.component(j) for j in range(3))
        for i in range(3)
    )


def extract_affine_poisson_geometry(
    graph: Graph, cell: basix.CellType
) -> tuple[Graph, GeometryKernelSpec | None]:
    """Extract the affine tetrahedral Poisson metric contraction from ``graph``.

    Unsupported expressions are deliberately left unchanged. This keeps the
    experimental path local to the one form and geometry mapping whose
    semantics are currently explicit.
    """
    if cell != basix.CellType.tetrahedron:
        return graph, None

    replacements: dict[GraphNode, GraphNode] = {}
    for node in graph.ordered_nodes():
        if (contraction := _poisson_metric_contraction(node)) is not None:
            replacements[node] = contraction

    if not replacements:
        return graph, None
    return replace(graph, replacements), GeometryKernelSpec()
