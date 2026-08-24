"""Form algorithms."""

from uflx.expressions import expression_sum, AbstractExpression
from uflx.graphs import Graph, generate_graph, GraphNode
from uflx.graphs.algorithms import reconstruct_node
from uflx.operators import Inner


def expand_inner_products(graph: Graph) -> Graph:
    """Replace inner products with sums over products of components."""
    new_nodes: dict[GraphNode, GraphNode] = {}
    for node in graph.ordered_nodes():
        if isinstance(node, Inner):
            assert node.first.value_shape == node.second.value_shape
            first = new_nodes[node.first]
            second = new_nodes[node.second]
            assert isinstance(first, AbstractExpression)
            assert isinstance(second, AbstractExpression)
            match len(node.first.value_shape):
                case 0:
                    new_nodes[node] = first * second
                case 1:
                    new_nodes[node] = expression_sum(
                        first.component(i) * second.component(i)
                        for i in range(first.value_shape[0])
                    )
                case _:
                    raise NotImplementedError()
        else:
            new_nodes[node] = reconstruct_node(node, new_nodes)
    return generate_graph(new_nodes[graph.root])
