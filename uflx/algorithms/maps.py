"""Algorithms related to push forward and pull back maps."""

from typing import Protocol, runtime_checkable

from uflx.algorithms.reconstruct import reconstruct_node
from uflx.graphs import GraphNode, as_graph


@runtime_checkable
class PullBackToReference(Protocol):
    """Pull a node back to the reference cell."""

    def pull_back_to_reference(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the reference cell."""


def pull_back_to_reference(
    expression: GraphNode,
) -> GraphNode:
    """Pull terms in integrals back to reference values."""
    node_map: dict[GraphNode, GraphNode] = {}
    for node in as_graph(expression).ordered_nodes():
        if isinstance(node, PullBackToReference):
            node_map[node] = node.pull_back_to_reference(node_map)
        elif any(a in node_map for a in node.successors):
            node_map[node] = reconstruct_node(node, node_map)

    return node_map.get(expression, expression)
