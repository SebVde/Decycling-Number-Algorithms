import networkx as nx
import random
import itertools


def construct_H(G, F, nb):
    H = nx.subgraph(G, nb)
    for u, v in itertools.combinations(nb, 2):
        if len(set(nx.common_neighbors(G, u, v)) - F) > 0:
            H.add_edge(u, v)

    return H


def is_foldable(G, v):
    nb = set(nx.neighbors(G, v))
    for x, y, z in itertools.combinations(nb, 3):
        if not (G.has_edge(x, y) or G.has_edge(y, z) or G.has_edge(z, x)):
            # If no edges between x-y-z, there is an anti-triangle so not foldable
            return False

    return True


def get_anti_edges(G, v):
    nb = set(nx.neighbors(G, v))
    anti_edges = set()
    for x, y in itertools.combinations(nb, 2):
        if not G.has_edge(x, y):
            anti_edges.add((x, y))

    return anti_edges


def fold_graph(G, v, anti_edges):
    nb_v = nx.neighbors(G, v)
    folded_G = nx.subgraph(G, G.nodes - nb_v - {v})
    added_nodes = set()

    for i, j in anti_edges:
        n = (i, j)
        folded_G.add_node(n)
        added_nodes.add(n)
        for u in set(nx.neighbors(G, i)) | set(nx.neighbors(G, j)):
            if u != v and u not in nb_v:
                folded_G.add_edge(n, u)

    for u, v in itertools.combinations(added_nodes, 2):
        folded_G.add_edge(u, v)

    return folded_G


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

    for u, v in itertools.combinations(G.nodes, 2):
        nb_u = set(nx.neighbors(G, u))
        nb_v = set(nx.neighbors(G, v))
        if (nb_v | {v}).issubset(nb_u | {u}):
            return get_max_indep_set(nx.subgraph(G, G.nodes - {u}))

    node_degrees = nx.degree(G)
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


def preprocess(G, F):
    if nx.number_connected_components(G) > 1:
        res = 0
        for c in nx.connected_components(G):
            if nx.is_forest(c):
                res += len(c.nodes)
            else:
                res += get_mif_len(c, set())

        return res

    sg_F = nx.subgraph(G, F)
    if len(sg_F.edges) != 0:
        pass


def get_mif_len(G, F):
    if F == set(G.nodes):
        return len(G.nodes)

    if len(F) == 0:
        max_degree = max(dict(nx.degree(G)).values())

        if max_degree < 2:
            return len(G.nodes)
        else:
            t = next(n for n, d in G.degree() if d >= 2)
            new_G = nx.subgraph(G, G.nodes - {t})
            return max(get_mif_len(G, F | {t}), get_mif_len(new_G, F))

    t = random.choice(list(F))
    nb = nx.neighbors(G, t)
    if G.nodes - F == nb:
        return len(F) + get_max_indep_set(construct_H(G, F - {t}, nb))

    pass


def get_decycling_number_mif_v2(G):
    # TODO preprocessing must be done at each call of get_mif_len
    if nx.is_forest(G):
        return 0

    if nx.number_connected_components(G) > 1:
        res = 0
        for c in nx.connected_components(G):
            if nx.is_forest(c):
                res += len(c.nodes)
            else:
                res += get_mif_len(c, set())

        return len(G.nodes) - res

    else:
        return len(G.nodes) - get_mif_len(G, set())
