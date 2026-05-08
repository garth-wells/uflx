"""Graphs."""

from __future__ import annotations

from typing import Protocol, Self
from networkx import DiGraph as Graph


class GraphNode(Protocol):
    """A node in a graph."""

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    def reconstruct(self, replacements: dict[GraphNode, GraphNode]) -> Self:
        """Reconstruct this node with some arguments replaced."""


class RepresentedByGraph(Protocol):
    """An object whose construction is represented by a graph."""

    @property
    def graph(self) -> Graph:
        """The graph that represents this object."""


def generate_graph(node: GraphNode) -> Graph:
    """Generate the graph that represents the construction of a node."""
    graph = Graph()

    added_nodes = {node}
    graph.add_node(node)
    while len(added_nodes) > 0:
        for n in added_nodes:
            for successor in n.successors:
                graph.add_node(successor)
                graph.add_edge(n, successor)
        added_nodes = set().union(*[n.successors for n in added_nodes])

    return graph
