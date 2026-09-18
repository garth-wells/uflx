"""Simplifying expressions."""

from uflx.expressions import AbstractExpression, Product
from typing import runtime_checkable, Protocol

from uflx.graphs import GraphNode, as_graph
from uflx.algorithms.reconstruct import reconstruct_node



@runtime_checkable
class SimplifiableInProduct(Protocol):
    """An expression that can be combined with others within a multiplication."""

    def simplified_product(self, other: AbstractExpression) -> AbstractExpression | None:
        """Return a single expression representing the simplified product.

        This function should return None if no simplification can be made.
        """


def simplify_product(items: Sequence[AbstractExpression]) -> AbstractExpression:
    """Simplify a product."""
    size = -1
    while len(items) != size:
        print(items, size)
        size = len(items)
        for i, item in enumerate(items):
            if isinstance(item, SimplifiableInProduct):
                for j, item2 in enumerate(items):
                    if i != j and (s := item.simplified_product(item2)) is not None:
                        items = [it for k, it in enumerate(items) if k not in [i, j]] + [s]
                        break
                else:
                    continue
                break

    if len(items) == 1:
        return items[0]
    else:
        return Product(items)


def simplify(expression: GraphNode) -> GraphNode:
    """Apply simplifications to an expression."""
    graph = as_graph(expression)
    assert graph.is_dag()

    node_map: dict[GraphNode, GraphNode] = {}
    for node in graph.ordered_nodes():
        if isinstance(node, Product):
            node_map[node] = simplify_product([node_map.get(i, i) for i in node._items])
        elif any(a in node_map for a in node.successors):
            node_map[node] = reconstruct_node(node, node_map)

    return node_map.get(graph.root, graph.root)
