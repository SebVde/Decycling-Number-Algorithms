import networkx as nx


def find_min_fvs(G):
    """
    Finds an approximate minimal FVS for the input graph.

    Args:
        G (nx.Graph): The input graph.

    Returns:
        set: A set of vertices forming the approximate minimal FVS.
    """
    fvs = set()
    while True:
        try:
            cycle = nx.find_cycle(G)
        except nx.NetworkXNoCycle:
            break

        cycle_nodes = {node for e in cycle for node in e}
        # Select the node with the highest degree to remove
        v = max(cycle_nodes, key=lambda x: G.degree(x))
        fvs.add(v)
        G.remove_node(v)

    return fvs


def approx_decycling_number_stanojevic(G):
    """
    Approximates the decycling number of the graph using Stanojevic's algorithm.

    Args:
        G (nx.Graph): The input graph.

    Returns:
        int: The approximate decycling number of the graph.
    """
    return len(find_min_fvs(G.copy()))
