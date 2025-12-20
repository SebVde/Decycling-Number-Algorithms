import networkx as nx
import itertools


def sort_nodes(node):
    """
    Function to sort nodes for consistent ordering.

    Args:
        node: A node in the graph, which can be of any hashable type.

    Returns:
        tuple: A tuple where the first element indicates if the node is a tuple,
               and the second element is the string representation of the node.
    """
    is_tuple = isinstance(node, tuple)
    string_rep = str(node)
    # False comes before True when sorting so non-tuple nodes are prioritized
    # If two nodes are of the same type, sort by their string representation
    return (is_tuple, string_rep)


def get_non_trivial_components(G):
    """
    Retrieves the connected components of G that contain more than one node.

    Args:
        G (networkx.Graph): A networkx graph.

    Returns:
        list: A list of sets reprenting the non-trivial components of G.
    """
    components = []
    for c in nx.connected_components(G):
        if len(c) > 1:
            components.append(c)

    return components


def construct_H(G, F, nb):
    """
    Constructs a subgraph H from the neighbors of a given set of nodes and adding edges
    between nodes that share a common neighbor in F.

    Args:
        G (networkx.Graph): A networkx graph.
        F (set): A set of nodes in G.
        nb (set): A set of neighbors in G.

    Returns:
        networkx.Graph: The constructed subgraph H.
    """
    H = nx.subgraph(G, nb).copy()
    for u, v in itertools.combinations(nb, 2):
        if not set(nx.common_neighbors(G, u, v)).isdisjoint(F):
            H.add_edge(u, v)

    return H


def is_foldable(G, v):
    """
    Checks if a node v in the graph is foldable. A node is foldable if none of its neighbors
    form an anti-triangle. An anti-triangle is a set of three nodes with no edges between them.

    Args:
        G (networkx.Graph): A networkx graph.
        v: A node in G.

    Returns:
        bool: True if the node is foldable, False otherwise.
    """
    nb = sorted(nx.neighbors(G, v), key=sort_nodes)
    for x, y, z in itertools.combinations(nb, 3):
        if not (G.has_edge(x, y) or G.has_edge(y, z) or G.has_edge(z, x)):
            # If no edges between x-y-z, there is an anti-triangle so not foldable
            return False

    return True


def get_anti_edges(G, v):
    """
    Finds all anti-edges among the neighbors of a node v. An anti-edge is a pair of
    nodes that are not connected by an edge.

    Args:
        G (networkx.Graph): A networkx graph.
        v: A node in G.

    Returns:
        set: A set of tuples representing anti-edges.
    """
    nb = nx.neighbors(G, v)
    anti_edges = set()
    for x, y in itertools.combinations(nb, 2):
        if not G.has_edge(x, y):
            anti_edges.add((x, y))

    return anti_edges


def fold_graph(G, v, anti_edges):
    """
    Folds the graph by contracting anti-edges of a node v and removing its neighbors.

    Args:
        G (networkx.Graph): A networkx graph.
        v: A node in G.
        anti_edges: A set of anti-edges to be contracted.

    Returns:
        networkx.Graph: The folded graph.
    """
    nb_v = set(nx.neighbors(G, v))
    folded_G = G.copy()
    added_nodes = set()

    for i, j in anti_edges:
        # Create a new node representing the contraction of i and j
        n = (i, j)
        folded_G.add_node(n)
        added_nodes.add(n)
        for u in set(nx.neighbors(G, i)) | set(nx.neighbors(G, j)):
            # Add edges between the new node and the neighbors of i and j except v and those in nb_v
            # (these nodes will be removed later)
            if u != v and u not in nb_v:
                folded_G.add_edge(n, u)

    for a, b in itertools.combinations(added_nodes, 2):
        # Add edges between each newly created node
        folded_G.add_edge(a, b)

    return nx.subgraph(folded_G, set(folded_G.nodes) - nb_v - {v})


def get_2_hop_neighbors(G, v):
    """
    Finds all nodes that are two hops away from a given node v.

    Args:
        G (networkx.Graph): A networkx graph.
        v: A node in G.

    Returns:
        set: A set of nodes that are two hops away from v.
    """
    first_neighbors = set(G.neighbors(v))
    second_neighbors = set()
    for u in first_neighbors:
        second_neighbors.update(G.neighbors(u))

    second_neighbors -= first_neighbors
    second_neighbors.discard(v)
    return second_neighbors


def is_complete(G):
    """
    Checks if a graph is complete. A graph is complete if every pair of nodes is connected.

    Args:
        G (networkx.Graph): A networkx graph.

    Returns:
        bool: True if the graph is complete, False otherwise.
    """
    return len(G.edges) == len(G.nodes) * (len(G.nodes) - 1) / 2


def get_mirrors(G, v):
    """
    Finds all mirror nodes of a given node v. A node u is a mirror of v if their
    neighborhoods are identical or if the subgraph induced by the difference of their
    neighborhoods is complete.

    Args:
        G (networkx.Graph): A networkx graph.
        v: A node in G.

    Returns:
        set: A set of mirror nodes of v.
    """
    mirrors = set()

    for u in get_2_hop_neighbors(G, v):
        nb_u = set(nx.neighbors(G, u))
        nb_v = set(nx.neighbors(G, v))
        if (len(nb_v - nb_u) == 0) or is_complete(nx.subgraph(G, nb_v - nb_u)):
            mirrors.add(u)

    return mirrors


def get_max_indep_set(G):
    """
    Computes the size of the maximum independent set of a graph.

    Args:
        G (networkx.Graph): A networkx graph.

    Returns:
        int: The size of the maximum independent set.
    """
    if len(G.nodes) == 0:
        return 0

    if nx.number_connected_components(G) > 1:
        res = 0
        for c in nx.connected_components(G):
            res += get_max_indep_set(G.subgraph(c))

        return res

    for u, v in itertools.combinations(sorted(G.nodes, key=sort_nodes), 2):
        nb_u = set(nx.neighbors(G, u))
        nb_v = set(nx.neighbors(G, v))
        if (nb_v | {v}).issubset(nb_u | {u}):
            return get_max_indep_set(nx.subgraph(G, set(G.nodes) - {u}))

    node_degrees = list(nx.degree(G))
    max_deg_4 = [(n, d) for n, d in node_degrees if d <= 4]
    sorted_deg = sorted(max_deg_4, key=lambda x: x[1])

    for comb in sorted_deg:
        v = comb[0]
        deg = comb[1]
        anti_edges = get_anti_edges(G, v)
        if (deg < 4 and is_foldable(G, v)) or (deg == 4 and len(anti_edges) < 4):
            return 1 + get_max_indep_set(fold_graph(G, v, anti_edges))

    v = max(node_degrees, key=lambda x: x[1])[0]
    return max(
        get_max_indep_set(nx.subgraph(G, G.nodes - {v} - get_mirrors(G, v))),
        1 + get_max_indep_set(nx.subgraph(G, G.nodes - {v} - set(nx.neighbors(G, v)))),
    )


def get_generalized_neighbors(G, F, active_v, v):
    """
    Finds the generalized neighbors of a node v considering a set of nodes F
    and an active node.

    Args:
        G (networkx.Graph): A networkx graph.
        F: A set of fixed nodes in G.
        active_v: The active node in G.
        v: The node for which generalized neighbors are computed.

    Returns:
        set: A set of generalized neighbors of v.
    """
    K = {u for u in (F - {active_v}) if G.has_edge(u, v)}
    new_G = G.copy()

    for n in K:
        new_G = nx.contracted_nodes(new_G, v, n, self_loops=False)

    gen_nb = set(new_G.neighbors(v)) - {active_v}
    return gen_nb


def get_mif_len(G, F, active_v):
    """
    Computes the length of a maximum induced forest (MIF) of a graph. This MIF must include all nodes in F.

    Args:
        G (networkx.Graph): A networkx graph.
        F: A set of fixed nodes in G.
        active_v: The active node in G.

    Returns:
        int: The length of a MIF of G.
    """
    if nx.number_connected_components(G) > 1:
        # If G is disconnected, compute MIF length for each component separately
        res = 0
        for c in nx.connected_components(G):
            if nx.is_forest(nx.subgraph(G, c)):
                res += len(c)
            else:
                res += get_mif_len(
                    nx.subgraph(G, c), F & set(c), active_v if active_v in c else None
                )

        return res

    sg_F = nx.subgraph(G, F)
    # Verify if F is acyclic
    if len(sg_F.nodes) > 0 and not nx.is_forest(sg_F):
        return 0

    if (
        len(sg_F.edges) != 0
    ):  # If F is not independent (if not every component of G[F] is an isolated vertex)

        new_G = G.copy()
        new_F = set(F)
        for T in get_non_trivial_components(sg_F):
            # Get all neighbors of T in G and need to remove those with more than 1 connection to T
            nb_T = set()
            for v in T:
                nb_T.update(set(G.neighbors(v)))
            nb_T -= T
            vertices_to_remove = set()
            for v in nb_T:
                connections = sum(1 for u in T if G.has_edge(u, v))
                if connections > 1:
                    vertices_to_remove.add(v)

            v_T = list(sorted(T, key=sort_nodes))[0]

            # Contract all nodes in T into v_T
            for n in T - {v_T}:
                new_G = nx.contracted_nodes(new_G, v_T, n, self_loops=False)

            new_G = nx.subgraph(new_G, new_G.nodes - vertices_to_remove)

            if (active_v is not None) and (active_v in T):
                active_v = v_T

            new_F -= T
            new_F.add(v_T)

        return get_mif_len(new_G, new_F, active_v) + len(F - new_F)

    else:
        return main_procedure(G, F, active_v)


def main_procedure(G, F, active_v):
    """
    Main procedure for computing the length of a maximum induced forest (MIF) of G.
    This MIF must include all nodes in F.

    Args:
        G (networkx.Graph): A networkx graph.
        F: A set of fixed nodes in G.
        active_v: The active node in G.

    Returns:
        int: The length of a MIF in G.
    """
    if F == set(G.nodes):
        return len(G.nodes)

    if len(F) == 0:
        max_degree = int(max(dict(nx.degree(G)).values()))

        if max_degree < 2:
            return len(G.nodes)
        else:
            # Choose a node t with degree at least 2
            # t is either contained in a MIF or not, thus we explore both subproblems
            t = list(sorted([n for n, d in G.degree() if d >= 2], key=sort_nodes))[0]
            new_G = nx.subgraph(G, G.nodes - {t})
            return max(
                get_mif_len(G, F | {t}, active_v), get_mif_len(new_G, F, active_v)
            )

    if active_v is None:
        active_v = list(sorted(F, key=sort_nodes))[0]

    nb = set(nx.neighbors(G, active_v))
    if set(G.nodes) - F == nb:
        # Get a maximum independent set of the constructed graph H
        return len(F) + get_max_indep_set(construct_H(G, F - {active_v}, nb))

    for v in sorted(nb, key=sort_nodes):
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        if len(gen_nb) < 2:
            # Include v in the MIF
            return get_mif_len(G, F | {v}, active_v)

    for v in sorted(nb, key=sort_nodes):
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        if len(gen_nb) > 3:
            # Explore both subproblems of including and excluding v from the MIF
            return max(
                get_mif_len(G, F | {v}, active_v),
                get_mif_len(nx.subgraph(G, G.nodes - {v}), F, active_v),
            )

    for v in sorted(nb, key=sort_nodes):
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        if len(gen_nb) == 2:
            # Either include v in the MIF or exclude v and include its two generalized neighbors in the MIF
            return max(
                get_mif_len(G, F | {v}, active_v),
                get_mif_len(
                    nx.subgraph(G, set(G.nodes) - {v}),
                    F | gen_nb,
                    active_v,
                ),
            )

    # If every v in nb has exactly 3 generalized neighbors, find a v that has at least
    # one generalized neighbor outside nb
    v = None
    good_gen_nb = None
    for n in nb:
        gen_nb = get_generalized_neighbors(G, F, active_v, n)
        if any(u not in nb for u in gen_nb):
            v = n
            good_gen_nb = gen_nb
            break

    if v is not None and good_gen_nb is not None:
        w1, w2, w3 = None, None, None
        for u in good_gen_nb:
            if w1 is None and u not in nb:
                w1 = u
            elif w2 is None:
                w2 = u
            elif w3 is None:
                w3 = u

            else:
                print("Problem")

        # Either include v in the MIF, or exclude v and include w1 (which is not in the neighborhood of
        # the active vertex), or exclude v and w1 and include w2 and w3 in the MIF
        return max(
            get_mif_len(G, F | {v}, active_v),
            get_mif_len(
                nx.subgraph(G, set(G.nodes) - {v}),
                F | {w1},
                active_v,
            ),
            get_mif_len(
                nx.subgraph(G, set(G.nodes) - {v, w1}),
                F | {w2, w3},
                active_v,
            ),
        )

    else:
        # Should not happen but for safety
        print("Problem")
        return Exception("No suitable v found")


def get_decycling_number_fomin(G):
    """
    Computes the decycling number of a graph using Fomin's algorithm.

    Args:
        G (networkx.Graph): A networkx graph.

    Returns:
        int: The decycling number of the graph.
    """
    if nx.is_forest(G):
        return 0
    return len(G.nodes) - get_mif_len(G, set(), None)
