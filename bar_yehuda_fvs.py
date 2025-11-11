import networkx as nx


def find_maximal_2_3_subgraph(G):
    pass


def subG_2_3(G):
    if nx.is_forest(G):
        return set()

    pass


def get_decycling_number_bar_yehuda(G):
    return len(subG_2_3(G))
