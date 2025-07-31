import networkx as nx


class CycleDetected(Exception):
    def __init__(self, message):
        super().__init__(message)


# Algorithm taken from the article "Exact Computation of Maximum Induced Forest"


def get_cnf(G, T):
    cnf = set()
    components = list(nx.connected_components(nx.subgraph(G, T)))

    for v in G.nodes():
        if v in T:
            continue

        # For each connected component from T, check if v has 2 or more neighbors in the component
        for comp in components:
            neighbors_in_comp = set(nx.neighbors(G, v)) & comp # & = intersection
            if len(neighbors_in_comp) >= 2:
                cnf.add(v)
                break

    return cnf


def T_update(G, T, K, v):
    pass


def K_update(G, T, K, W):
    pass


def find_mif(G, T, K):
    pass


def main_mif(G, K):
    K_induced = nx.subgraph(G, K)

    if nx.is_forest(K_induced):
        sg_nodes = G.nodes - get_cnf(G, K)
        sg = nx.subgraph(G, sg_nodes)
        find_mif(sg, set(), K)

    else:
        raise CycleDetected("The subgraph induced by the given K should NOT contain cycles")

