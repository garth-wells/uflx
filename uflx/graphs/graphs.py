"""Graphs."""

from __future__ import annotations

from typing import Protocol, Self
from networkx import DiGraph


class Graph(DiGraph):
    """An acyclic directed graph."""

    _root: GraphNode | None

    def __init__(self, *args, **kwargs):
        """Initialise."""
        self._root = None
        super().__init__(*args, **kwargs)

    @property
    def root(self) -> GraphNode:
        """Get the root node of the graph."""
        assert self._root is not None
        return self._root

    def set_root(self, node: GraphNode):
        """Set the root node of the graph."""
        self._root = node

    def add_root_node(self, node: GraphNode, *args, **kwargs):
        """Add a new node to the graph and set it to be the root."""
        self.add_node(node, *args, **kwargs)
        self._root = node


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
    graph.add_root_node(node)
    while len(added_nodes) > 0:
        for n in added_nodes:
            for successor in n.successors:
                graph.add_node(successor)
                graph.add_edge(n, successor)
        added_nodes = set().union(*[n.successors for n in added_nodes])

    return graph
