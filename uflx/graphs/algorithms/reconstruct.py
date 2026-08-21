"""Reconstruct."""
from typing import Any

from uflx.graphs.graphs import GraphNode


def apply_replacements(arg: Any, replacements: dict[GraphNode, GraphNode]) -> Any:
    """Apply replacements to an init arg."""
    if isinstance(arg, GraphNode):
        return replacements.get(arg, arg)
    if isinstance(arg, list):
        return [apply_replacements(i, replacements) for i in arg]
    if isinstance(arg, tuple):
        return tuple(apply_replacements(i, replacements) for i in arg)
    if isinstance(arg, dict):
        return {
            apply_replacements(i, replacements): apply_replacements(j, replacements)
            for i, j in arg.items()
        }

    return arg



def reconstruct_node(node: GraphNode, replacements: dict[GraphNode, GraphNode]) -> GraphNode:
    """Reconstruct a node with replacements made."""
    args = [apply_replacements(i, replacements) for i in node.init_args]
    return node.__class__(*args)
