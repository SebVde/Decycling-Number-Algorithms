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


def main_procedure(G, F, active_v):
    pass


def get_mif_len(G, F, active_v):
    if not nx.subgraph(G, F).is_forest():
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
                if new_G.number_of_edges(u, n) == 2:
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
