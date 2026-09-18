"""Simplifying expressions."""

from collections.abc import Sequence
from typing import Protocol, runtime_checkable

from uflx.algorithms.reconstruct import reconstruct_node
from uflx.expressions import AbstractExpression, Product, Sum
from uflx.graphs import GraphNode, as_graph


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


@runtime_checkable
class SimplifiableInSum(Protocol):
    """An expression that can be combined with others within a multiplication."""

    def simplified_sum(self, other: AbstractExpression) -> AbstractExpression | None:
        """Return a single expression representing the simplified sum.

        This function should return None if no simplification can be made.
        """


def simplify_sum(items: Sequence[AbstractExpression]) -> AbstractExpression:
    """Simplify a sum."""
    size = -1
    while len(items) != size:
        size = len(items)
        for i, item in enumerate(items):
            if isinstance(item, SimplifiableInSum):
                for j, item2 in enumerate(items):
                    if i != j and (s := item.simplified_sum(item2)) is not None:
                        items = [it for k, it in enumerate(items) if k not in [i, j]] + [s]
                        break
                else:
                    continue
                break

    if len(items) == 1:
        return items[0]
    else:
        return Sum(items)


def _mapped_items(
    items: tuple[AbstractExpression, ...], node_map: dict[GraphNode, GraphNode]
) -> list[AbstractExpression]:
    """Get list of mapped items."""
    out = []
    for i in items:
        j = node_map.get(i, i)
        assert isinstance(j, AbstractExpression)
        out.append(j)
    return out


def simplify(expression: GraphNode) -> GraphNode:
    """Apply simplifications to an expression."""
    graph = as_graph(expression)
    assert graph.is_dag()

    node_map: dict[GraphNode, GraphNode] = {}
    for node in graph.ordered_nodes():
        if isinstance(node, Product):
            node_map[node] = simplify_product(_mapped_items(node._items, node_map))
        elif isinstance(node, Sum):
            node_map[node] = simplify_sum(_mapped_items(node._items, node_map))
        elif any(a in node_map for a in node.successors):
            node_map[node] = reconstruct_node(node, node_map)

    return node_map.get(graph.root, graph.root)
