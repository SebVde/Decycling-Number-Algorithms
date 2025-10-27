import networkx as nx


def get_mif_len(G, F):
    pass


def get_decycling_number_mif_v3(G):
    if nx.is_forest(G):
        return 0

    return len(G.nodes) - get_mif_len(G, set())
