import networkx as nx
from pyvis.network import Network

nt = nx.complete_graph(["A", "B", "C", "D"])

graph = Network()
graph.inherit_edge_colors(False)
graph.from_nx(nt)
graph.get_node("A")["color"] = "green"

graph.show("graph.html", notebook=False)