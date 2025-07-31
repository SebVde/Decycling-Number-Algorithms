import networkx as nx
from pyvis.network import Network
from exact_mif import main_mif

nt = nx.Graph()
nt.add_nodes_from(["A", "B", "C", "D"])
nt.add_edges_from([("A", "B"), ("B", "C"), ("C", "D")])

main_mif(nt, nt.nodes())

# graph = Network()
# graph.inherit_edge_colors(False)
# graph.from_nx(nt)
# graph.get_node("A")["color"] = "green"
#
# graph.show("graph.html", notebook=False)