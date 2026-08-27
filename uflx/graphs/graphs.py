"""Graphs."""

from __future__ import annotations

from collections.abc import Hashable, Iterable
from enum import Enum
from typing import Any, Protocol, cast, runtime_checkable

import networkx as nx


class NodeOrder(Enum):
    """A node order."""

    roots_first = 0
    leaves_first = 1


def _print_tree(graph: Graph, node: GraphNode, prefix: str, is_last: bool):
    connector = "└── " if is_last else "├── "
    print(prefix + connector + repr(node))
    successors = list(graph.successors(node))
    extension = "    " if is_last else "│   "
    for i, child in enumerate(successors):
        _print_tree(graph, child, prefix + extension, i == len(successors) - 1)


def print_node(graph: Graph, node: GraphNode, indentation: int = 0):
    """Print a graph using the node as the root node."""
    indent = "    " * indentation
    print(indent + repr(node))
    successors = list(graph.successors(node))
    for i, child in enumerate(successors):
        _print_tree(graph, child, indent, i == len(successors) - 1)


@runtime_checkable
class GraphNode(Hashable, Protocol):
    """A node in a graph."""

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""


class Graph(nx.DiGraph):
    """An acyclic directed graph."""

    _root: GraphNode | None

    def __init__(self, *args, **kwargs):
        """Initialise."""
        self._root = None
        super().__init__(*args, **kwargs)
        assert nx.is_directed_acyclic_graph(self)

    @property
    def root(self) -> GraphNode:
        """Get the root node of the graph."""
        assert self._root is not None
        return self._root

    @property
    def has_root(self) -> bool:
        """Check if this graph has a root node."""
        return self._root is not None

    def set_root(self, node: GraphNode):
        """Set the root node of the graph."""
        self._root = node

    def add_root_node(self, node: GraphNode, *args, **kwargs):
        """Add a new node to the graph and set it to be the root."""
        self.add_node(node, *args, **kwargs)
        self._root = node

    def print(self):
        """Print a graph."""
        print_node(self, self.root)

    def subgraph(self, node: GraphNode) -> Graph:  # type: ignore[override]
        """Get the subgraph with the input node as the root node."""
        return generate_graph(node)

    def ordered_nodes(self, order: NodeOrder = NodeOrder.leaves_first) -> Iterable[GraphNode]:
        """Iterate through the ordered graph nodes."""
        match order:
            case NodeOrder.roots_first:
                return cast(Iterable[GraphNode], nx.topological_sort(self))
            case NodeOrder.leaves_first:
                return cast(Iterable[GraphNode], reversed(list(nx.topological_sort(self))))
            case _:
                raise ValueError("Invalid node order")

    def descendants(self, node: GraphNode) -> set[GraphNode]:
        """Get all descendants of a node."""
        return cast(set[GraphNode], nx.descendants(self, node))

    def is_dag(self) -> bool:
        """Check if this graph is a directed acyclic graph."""
        return nx.is_directed_acyclic_graph(self)


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
