import networkx as nx


def find_min_fvs(G):
    fvs = set()
    while True:
        try:
            cycle = nx.find_cycle(G)
        except nx.NetworkXNoCycle:
            break

        cycle_nodes = {node for e in cycle for node in e}
        v = max(cycle_nodes, key=lambda x: G.degree(x))
        fvs.add(v)
        G.remove_node(v)

    return fvs


def approx_decycling_number_stanojevic(G):
    return len(find_min_fvs(G.copy()))
