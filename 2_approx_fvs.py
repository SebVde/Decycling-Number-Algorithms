import networkx as nx


def cleanup(G):
    G = G.copy()
    for node in list(G.nodes):
        if G.degree(node) <= 1:
            G.remove_node(node)
    return G

def get_fvs(G):
    F = set()

    for node in G.nodes:
        if G.nodes[node].get("weight", 1) == 0:
            F.add(node)

    G.remove_nodes_from(F)

    G = cleanup(G)
    i = 0

    while G.order() > 0:
        i += 1

        # a cycle C is semidisjoint if, for every vertex u of C, d(u) = 2 with at most one exception.
        # if G contains a semidisjoint cycle C then we can remove all vertices of C from G

        # C = nx.cycle_basis(G)
