import networkx as nx
from pyvis.network import Network
from exact_mif import get_mif_len
from naive import get_decycling_number

nt = nx.erdos_renyi_graph(12, 0.89)

# nt = nx.Graph()
# nt.add_nodes_from(["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"])
# nt.add_edges_from([
#     ("V1", "V2"),
#     ("V1", "V3"),
#     ("V2", "V4"),
#     ("V4", "V5"),
#     ("V3", "V5"),
#     ("V4", "V6"),
#     ("V6", "V7"),
#     ("V5", "V7"),
#     ("V7", "V8"),
#     ("V6", "V8"),
# ])

print(get_mif_len(nt, set()))
print(get_decycling_number(nt))

# graph = Network()
# graph.inherit_edge_colors(False)
# graph.options.edges.smooth.enabled = False
# graph.from_nx(nt)
# graph.get_node("V6")["color"] = "green"
#
# graph.show("graph.html", notebook=False)
