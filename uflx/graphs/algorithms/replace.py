"""Algorithm to replace node with different nodes."""

import networkx as nx

from uflx.graphs.graphs import Graph, GraphNode, generate_graph


def replace(graph: Graph, replacements: dict[GraphNode, GraphNode]) -> Graph:
    """Replace nodes in a graph.

    Args:
        graph: The graph
        replacements: A map from nodes to the nodes they should be replaced with

    Returns:
        A new graph with replacements made
    """
    assert nx.is_directed_acyclic_graph(graph)

    node_map: dict[GraphNode, GraphNode] = {}
    for node in reversed(list(nx.topological_sort(graph))):
        if node in replacements:
            node_map[node] = replacements[node]
        elif any(a in node_map for a in node.successors):
            node_map[node] = node.reconstruct(node_map)

    return generate_graph(node_map.get(graph.root, graph.root))
