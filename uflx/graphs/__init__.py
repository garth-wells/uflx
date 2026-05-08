"""Graphs."""

from networkx import DiGraph as Graph, union, is_directed_acyclic_graph as is_dag
from uflx.graphs.graphs import GraphNode, RepresentedByGraph, generate_graph
from uflx.graphs import algorithms
