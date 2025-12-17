import networkx as nx
import itertools


def sort_nodes(node):
    is_tuple = isinstance(node, tuple)
    safe_val = str(node)
    return (is_tuple, safe_val)


def get_non_trivial_components(G):
    components = []
    for c in nx.connected_components(G):
        if len(c) > 1:
            components.append(c)

    return components


def construct_H(G, F, nb):
    H = nx.subgraph(G, nb).copy()
    for u, v in itertools.combinations(nb, 2):
        if not set(nx.common_neighbors(G, u, v)).isdisjoint(F):
            H.add_edge(u, v)

    return H


def is_foldable(G, v):
    nb = sorted(nx.neighbors(G, v), key=sort_nodes)
    for x, y, z in itertools.combinations(nb, 3):
        if not (G.has_edge(x, y) or G.has_edge(y, z) or G.has_edge(z, x)):
            # If no edges between x-y-z, there is an anti-triangle so not foldable
            return False

    return True


def get_anti_edges(G, v):
    nb = nx.neighbors(G, v)
    anti_edges = set()
    for x, y in itertools.combinations(nb, 2):
        if not G.has_edge(x, y):
            anti_edges.add((x, y))

    return anti_edges


def fold_graph(G, v, anti_edges):
    nb_v = set(nx.neighbors(G, v))
    folded_G = G.copy()
    added_nodes = set()

    for i, j in anti_edges:
        n = (i, j)
        folded_G.add_node(n)
        added_nodes.add(n)
        for u in set(nx.neighbors(G, i)) | set(nx.neighbors(G, j)):
            if u != v and u not in nb_v:
                folded_G.add_edge(n, u)

    for a, b in itertools.combinations(added_nodes, 2):
        folded_G.add_edge(a, b)

    return nx.subgraph(folded_G, set(folded_G.nodes) - nb_v - {v})


def get_2_hop_neighbors(G, v):
    first_neighbors = set(G.neighbors(v))
    second_neighbors = set()
    for u in first_neighbors:
        second_neighbors.update(G.neighbors(u))

    second_neighbors -= first_neighbors
    second_neighbors.discard(v)
    return second_neighbors


def is_complete(G):
    return len(G.edges) == len(G.nodes) * (len(G.nodes) - 1) / 2


def get_mirrors(G, v):
    mirrors = set()

    for u in get_2_hop_neighbors(G, v):
        nb_u = set(nx.neighbors(G, u))
        nb_v = set(nx.neighbors(G, v))
        if (len(nb_v - nb_u) == 0) or is_complete(nx.subgraph(G, nb_v - nb_u)):
            mirrors.add(u)

    return mirrors


def get_max_indep_set(G):
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
    K = {u for u in (F - {active_v}) if G.has_edge(u, v)}
    new_G = G.copy()

    for n in K:
        new_G = nx.contracted_nodes(new_G, v, n, self_loops=False)

    gen_nb = set(new_G.neighbors(v)) - {active_v}
    return gen_nb


def get_mif_len(G, F, active_v):
    if nx.number_connected_components(G) > 1:
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
    # Verify is F is acyclic
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

    if F == set(G.nodes):
        return len(G.nodes)

    if len(F) == 0:
        max_degree = int(max(dict(nx.degree(G)).values()))

        if max_degree < 2:
            return len(G.nodes)
        else:
            t = list(sorted([n for n, d in G.degree() if d >= 2], key=sort_nodes))[0]
            new_G = nx.subgraph(G, G.nodes - {t})
            return max(
                get_mif_len(G, F | {t}, active_v), get_mif_len(new_G, F, active_v)
            )

    if active_v is None:
        active_v = list(sorted(F, key=sort_nodes))[0]

    nb = set(nx.neighbors(G, active_v))
    if set(G.nodes) - F == nb:
        return len(F) + get_max_indep_set(construct_H(G, F - {active_v}, nb))

    for v in sorted(nb, key=sort_nodes):
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        if len(gen_nb) < 2:
            return get_mif_len(G, F | {v}, active_v)

    for v in sorted(nb, key=sort_nodes):
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        if len(gen_nb) > 3:
            return max(
                get_mif_len(G, F | {v}, active_v),
                get_mif_len(nx.subgraph(G, G.nodes - {v}), F, active_v),
            )

    for v in sorted(nb, key=sort_nodes):
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        if len(gen_nb) == 2:
            return max(
                get_mif_len(G, F | {v}, active_v),
                get_mif_len(
                    nx.subgraph(G, set(G.nodes) - {v}),
                    F | gen_nb,
                    active_v,
                ),
            )

    # If every v in nb has exactly 3 generalized neighbors, find a v that has at least one generalized neighbor outside nb
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
        print("Problem")
        return Exception("No suitable v found")


def get_decycling_number_fomin(G):
    if nx.is_forest(G):
        return 0

    return len(G.nodes) - get_mif_len(G, set(), None)
