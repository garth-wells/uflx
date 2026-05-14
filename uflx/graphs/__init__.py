"""Graphs."""

from networkx import is_directed_acyclic_graph as is_dag
from networkx import union

from uflx.graphs import algorithms
from uflx.graphs.graphs import Graph, GraphNode, RepresentedByGraph, generate_graph
