import networkx as nx
from itertools import combinations


def naive_fvs(graph):
    all_nodes = set(graph.nodes)

    for i in range(1, len(all_nodes)):
        subsets = combinations(all_nodes, i)
        for subset in subsets:
            subset = set(subset)
            subgraph = nx.subgraph(graph, all_nodes - subset)
            if nx.is_forest(subgraph):
                return subset

    return None


def get_decycling_number(graph):
    return len(naive_fvs(graph))
