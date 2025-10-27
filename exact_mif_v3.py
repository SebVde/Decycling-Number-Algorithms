import itertools

import networkx as nx


def get_generalized_neighbors(G, F, active_v, v):
    gen_nb = set(nx.neighbors(G, v)) - F
    ens = set(nx.neighbors(G, v)) & (F - {active_v})
    for s in ens:
        gen_nb.update(set(nx.neighbors(G, s)) - F)

    return gen_nb


def get_non_trivial_components(G):
    components = []
    for c in nx.connected_components(G):
        if len(c) > 1:
            components.append(c)

    return components


def find_short_pair(G, F, active_v):
    for u, v in itertools.combinations(set(G.nodes()) - F, 2):
        nb_u = set(nx.neighbors(G, u))
        nb_v = set(nx.neighbors(G, v))
        # If parallel edges between u and v (so u,v is a short-cycle) or they have a common neighbor in F (also short-cycle)
        if (G.number_of_edges(u, v) > 1) or (
            G.number_of_edges(u, v) == 1 and any(w in F for w in nb_u & nb_v)
        ):
            if (active_v is None) or (active_v not in nb_u & nb_v):
                return u, v

    return None, None


def main_procedure(G, F, active_v):
    a, b = find_short_pair(G, F, active_v)
    if a is not None and b is not None:
        return max(
            get_mif_len(nx.subgraph(G, set(G.nodes) - {a}), F | {b}, active_v),
            get_mif_len(nx.subgraph(G, set(G.nodes) - {b}), F, active_v),
        )

    cut_v = set(nx.articulation_points(G))
    if len(cut_v) > 0:
        v = next(iter(cut_v))
        components = set(nx.connected_components(nx.subgraph(G, set(G.nodes) - {v})))
        H = min(components, key=len)

        F1 = F & (set(H) | {v})
        G1 = nx.subgraph(G, set(H) | {v})
        F2 = F - set(H)
        G2 = nx.subgraph(G, set(G.nodes) - set(H))
        F1_star = F1 | {v}

        if v in F:
            return get_mif_len(
                G1, F1, active_v if active_v in F1 else None
            ) + get_mif_len(G2, F2, active_v if active_v in F2 else None)

        sg_H = nx.subgraph(G, H)
        S1 = get_mif_len(sg_H, F & set(H), active_v if active_v in F & set(H) else None)
        S1_star = get_mif_len(G1, F1_star, active_v if active_v in F1_star else None)

        if S1_star > S1:
            return (
                S1_star - 1 + get_mif_len(G2, F2, active_v if active_v in F2 else None)
            )

        else:
            return get_mif_len(
                sg_H, F1, active_v if active_v in F1 else None
            ) + get_mif_len(
                nx.subgraph(G, set(G.nodes) - (set(H) | {v})),
                F - (set(H) | {v}),
                active_v if active_v in F - (set(H) | {v}) else None,
            )

    if len(F) == 0:
        degrees = dict(nx.degree(G))
        v = max(degrees, key=degrees.get)
        new_G = nx.subgraph(G, G.nodes - {v})
        return max(get_mif_len(G, F | {v}, v), get_mif_len(new_G, F, active_v))

    if active_v is None:
        active_v = next(iter(F))

    pass


def get_mif_len(G, F, active_v):
    if len(F) > 1 and not nx.subgraph(G, F).is_forest():
        print("Can't reduce cause F is not acyclic")
        return Exception

    new_G = G.copy()
    new_F = set(F)
    S = set()

    while True:
        # Step 1
        sg_F = nx.subgraph(new_G, new_F)
        non_trivial_components = get_non_trivial_components(sg_F)
        if len(non_trivial_components) > 0:
            T = non_trivial_components[0]
            v = next(iter(T))
            if (active_v is not None) and (active_v in T):
                v = active_v

            for n in T - {v}:
                new_G = nx.contracted_nodes(new_G, v, n, self_loops=False)

            S = S | (set(T) - {v})
            new_F = new_F - (set(T) - {v})
            continue

        # Step 2
        v = None
        # If we find a node v not in new_F that has 2 parallel edges to a node u in new_F, we remove it from new_G
        for n in set(new_G.nodes()) - new_F:
            nb_in_F = set(new_G.neighbors(n)) & new_F
            for u in nb_in_F:
                if new_G.number_of_edges(u, n) > 1:
                    v = n
                    break
            if v is not None:
                break

        if v is not None:
            new_G.remove_node(v)
            continue

        # Step 3
        for n, deg in new_G.degree():
            if deg == 1:
                v = n
                break

        if v is not None:
            S.add(v)
            new_G.remove(v)
            continue

        # Step 4
        for n, deg in new_G.degree():
            if n not in new_F and (
                deg == 2
                or len(get_generalized_neighbors(new_G, new_F, active_v, n)) <= 1
            ):
                v = n
                break

        if v is not None:
            new_F.add(v)
            continue

    if len(new_G.nodes) == 0:
        return len(S)
    else:
        return len(S) + main_procedure(new_G, new_F, active_v)


def get_decycling_number_mif_v3(G):
    if nx.is_forest(G):
        return 0

    return len(G.nodes) - get_mif_len(G, set(), None)
