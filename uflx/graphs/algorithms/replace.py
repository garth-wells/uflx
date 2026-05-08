"""Algorithm to replace node with different nodes."""

import networkx as nx

from uflx.graphs.graphs import Graph, GraphNode


def replace(graph: Graph, replacements: dict[GraphNode, GraphNode]) -> Graph:
    """Replace nodes in a graph.

    Args:
        graph: The graph
        replacements: A map from nodes to the nodes they should be replaced with

    Returns:
        A new graph with replacements made
    """
    assert nx.is_directed_acyclic_graph(graph)

    new_graph = Graph()

    node_map: dict[GraphNode, GraphNode] = {}
    for node in reversed(list(nx.topological_sort(graph))):
        if node in replacements:
            new_node = replacements[node]
        else:
            new_node = node.reconstruct(node_map)
        if new_node != node:
            node_map[node] = new_node
        new_graph.add_node(new_node)
        for s in graph.successors(node):
            new_graph.add_edge(new_node, node_map.get(s, s))

    new_graph.set_root(node_map.get(graph.root, graph.root))

    return new_graph
