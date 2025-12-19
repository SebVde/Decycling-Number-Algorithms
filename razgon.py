import networkx as nx


def get_cnf(G, T):
    """
    Computes the set of conflicting vertices in G that are not in T but have at least 2 neighbors
    in the same connected component of the subgraph of G induced by T.

    Args:
        G (networkx.Graph): A networkx graph.
        T (set): A subset of vertices in G.

    Returns:
        set: A set of vertices representing the conflicting vertices.
    """
    cnf = set()
    components = list(nx.connected_components(nx.subgraph(G, T)))

    if len(components) > 0:
        for v in sorted(set(G.nodes) - T):
            # For each connected component from T, check if v has 2 or more neighbors in the component
            for comp in components:
                neighbors_in_comp = set(nx.neighbors(G, v)) & comp
                if len(neighbors_in_comp) >= 2:
                    cnf.add(v)
                    break  # v only needs 2 neighbors in one component to be conflicting

    return cnf


def get_bnd(G, T):
    """
    Computes the boundary vertices of T in G, i.e., the vertices that are not in T but have exactly
    one neighbor in T.

    Args:
        G (networkx.Graph): A networkx graph.
        T (set): A subset of vertices in G.

    Returns:
        set: A set of vertices representing the boundary vertices of T in G.
    """
    bnd = set()

    for v in set(G.nodes) - T:
        neighbors_to_T = set(nx.neighbors(G, v)) & T
        if len(neighbors_to_T) == 1:
            bnd.add(v)

    return bnd


def get_K_connected(G, K, v):
    """
    Computes the subset of vertices that are K-connected to v. A vertex u is K-connected to v if
    it is adjacent to v or has a path to v such that the whole path is in K.

    Args:
        G (networkx.Graph): A networkx graph.
        K (set): A subset of vertices in G.
        v: A vertex in G.

    Returns:
        set: A set of vertices representing the K-connected vertices to v.
    """
    result = set(nx.neighbors(G, v))
    nb_in_K = result & K
    sg_K = nx.subgraph(G, K)
    if len(nb_in_K) > 0:
        for u in nb_in_K:
            result.update(nx.node_connected_component(sg_K, u))

    if v in K:
        result.update(nx.node_connected_component(sg_K, v))

    return result


def T_update(G, T, K, v):
    """
    Updates the sets T and K by adding the K-connected vertices to T and removing them from K,
    and removes conflicting vertices from the graph.

    Args:
        G (networkx.Graph): A networkx graph.
        T (set): A subset of vertices in G.
        K (set): A subset of vertices in G.
        v: A vertex in G.

    Returns:
        tuple: A tuple containing the updated graph, T, and K.
    """
    S = get_K_connected(G, K, v) & K
    new_T = T | S | {v}  # Update T with the K-connected vertices and v
    new_K = K - (S | {v})  # Update K by removing the K-connected vertices and v
    sg_nodes = set(G.nodes) - get_cnf(G, new_T)
    sg = nx.subgraph(
        G, sg_nodes
    )  # Create a subgraph induced by the nodes that are not conflicting with the new T
    return sg, new_T, new_K


def K_update(G, T, K, W):
    """
    Updates the graph G, T, and K based on the set W.

    Args:
        G (networkx.Graph): A networkx graph.
        T (set): A subset of vertices in G.
        K (set): A subset of vertices in G.
        W (set): A subset of 2 vertices in G.

    Returns:
        tuple or None: A tuple containing the updated graph, T, and K, or None if the subgraph
        induced by T, K, and W contains cycles.
    """
    sg = nx.subgraph(G, T | K | W)

    if nx.is_forest(sg):
        bnd = get_bnd(G, T)
        if len(W & bnd) == 0:
            sg_nodes = set(G.nodes) - get_cnf(G, K | W)
            return nx.subgraph(G, sg_nodes), T, K | W

        else:
            v1 = next(iter(W & bnd))
            new_G, new_T, new_K = T_update(G, T, K, v1)
            v2 = set(W - {v1}).pop()  # W has always size 2

            if v2 in get_bnd(new_G, new_T):
                return T_update(new_G, new_T, new_K, v2)
            else:
                sg_nodes = set(new_G.nodes) - get_cnf(new_G, new_K | {v2})
                return nx.subgraph(new_G, sg_nodes), new_T, new_K | {v2}

    else:
        return None  # If the subgraph contains cycles


def find_mif_len(G, T, K):
    """
    Searches for the length of a (T union K)-MIF of the graph G. A (T union K)-MIF is a maximal
    induced forest that includes all vertices of T and K.

    Args:
        G (networkx.Graph): A networkx graph.
        T (set): A subset of vertices in G.
        K (set): A subset of vertices in G.

    Returns:
        int: The length of a (T union K)-MIF of G.
    """
    if T | K == set(G.nodes):
        return len(T | K)

    else:
        v = None
        bnd = get_bnd(G, T)
        if len(bnd) > 0:
            v = list(sorted(bnd))[0]
        else:
            v = list(sorted(set(G.nodes) - (T | K)))[0]

        new_G, new_T, new_K = T_update(G, T, K, v)
        new_S = find_mif_len(new_G, new_T, new_K)

        W = get_K_connected(G, K, v) & set(G.nodes) - (T | K | {v})

        match len(W):
            case 0 | 1:
                return new_S
            case 2:
                sg = nx.subgraph(G, set(G.nodes) - {v})
                result = K_update(sg, T, K, W)
                if result is not None:
                    return max(new_S, find_mif_len(*result))
                else:
                    return new_S
            case _:
                return max(
                    new_S,
                    find_mif_len(nx.subgraph(G, set(G.nodes) - {v}), T, K),
                )


def get_decycling_number_razgon(G):
    """
    Computes the decycling number of a graph G using Razgon's algorithm.

    Args:
        G (networkx.Graph): A networkx graph.

    Returns:
        int: The decycling number of G.
    """
    if nx.is_forest(G):
        return 0
    return len(set(G.nodes)) - find_mif_len(G, set(), set())
