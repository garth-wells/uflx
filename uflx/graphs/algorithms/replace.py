"""Algorithm to replace node with different nodes."""

from uflx.graphs.algorithms.reconstruct import reconstruct_node
from uflx.graphs.graphs import Graph, GraphNode, generate_graph, NodeOrder


def replace(graph: Graph, replacements: dict[GraphNode, GraphNode]) -> Graph:
    """Replace nodes in a graph.

    Args:
        graph: The graph
        replacements: A map from nodes to the nodes they should be replaced with

    Returns:
        A new graph with replacements made
    """
    assert graph.is_dag()

    node_map: dict[GraphNode, GraphNode] = {}
    for node in graph.ordered_nodes():
        if node in replacements:
            node_map[node] = replacements[node]
        elif any(a in node_map for a in node.successors):
            node_map[node] = reconstruct_node(node, node_map)

    return generate_graph(node_map.get(graph.root, graph.root))
