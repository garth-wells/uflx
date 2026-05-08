"""Graphs."""

from networkx import DiGraph as Graph
from networkx import is_directed_acyclic_graph as is_dag
from networkx import union

from uflx.graphs import algorithms
from uflx.graphs.graphs import GraphNode, RepresentedByGraph, generate_graph
