import itertools
import networkx as nx


def get_generalized_neighbors(G, F, active_v, v):
    gen_nb = set(nx.neighbors(G, v)) - F
    excl = {active_v} if active_v is not None else set()
    ens = set(nx.neighbors(G, v)) & (F - excl)
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
    for u, v in itertools.combinations(set(G.nodes) - F, 2):
        nb_u = set(nx.neighbors(G, u))
        nb_v = set(nx.neighbors(G, v))
        # If parallel edges between u and v (so u,v is a short-cycle) or they have a common neighbor in F (also short-cycle)
        if (G.number_of_edges(u, v) > 1) or (
            G.number_of_edges(u, v) == 1 and any(w in F for w in nb_u & nb_v)
        ):
            if not ((active_v is not None) and (active_v in nb_u or active_v in nb_v)):
                return u, v

    return None, None


def is_trigger_vertex(G, F, active_v, v):
    nb_active = set(nx.neighbors(G, active_v))
    for u in nb_active - F:
        gen_nb = get_generalized_neighbors(G, F, active_v, u)
        if len(gen_nb - nb_active) >= 3 and v in gen_nb:
            nb_v = set(nx.neighbors(G, v))
            nb_u = set(nx.neighbors(G, u))
            s_set = F - nb_v
            v_prime_set = F & nb_v
            for s in s_set:
                if nb_u == {active_v, v, s}:
                    return True

                for v_prime in v_prime_set:
                    d_v_prime = G.degree(v_prime)
                    if d_v_prime == 2 and (
                        nb_u == {active_v, v_prime, s}
                        or nb_u == {active_v, v, v_prime, s}
                    ):
                        return True

    return False


def find_optimal_v(G, F, active_v):
    nb_active = set(nx.neighbors(G, active_v))
    G_not_F = set(G.nodes) - F
    possible = set()
    for v in G_not_F:
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        if len(gen_nb) < 3:
            continue

        if is_trigger_vertex(G, F, active_v, v):
            return v

        if len(gen_nb) == 3:
            possible.add(v)

    if len(possible) > 0:
        return max(
            possible,
            key=lambda x: len(get_generalized_neighbors(G, F, active_v, x) - nb_active),
        )

    else:
        return max(
            G_not_F,
            key=lambda x: (
                len(set(nx.neighbors(G, x)) & (F - {active_v})),
                len(get_generalized_neighbors(G, F, active_v, x) & nb_active),
                len(get_generalized_neighbors(G, F, active_v, x)),
            ),
        )


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
        components = list(nx.connected_components(nx.subgraph(G, set(G.nodes) - {v})))
        H = min(components, key=lambda x: len(x))

        F1 = F & (H | {v})
        G1 = nx.subgraph(G, H | {v})
        F2 = F - H
        G2 = nx.subgraph(G, set(G.nodes) - H)
        F1_star = F1 | {v}

        if v in F:
            return (
                get_mif_len(G1, F1, active_v if active_v in F1 else None)
                + get_mif_len(G2, F2, active_v if active_v in F2 else None)
                - 1
            )  # Because v counted twice

        sg_H = nx.subgraph(G, H)
        S1 = get_mif_len(sg_H, F & H, active_v if active_v in F & H else None)
        S1_star = get_mif_len(G1, F1_star, active_v if active_v in F1_star else None)

        if S1_star > S1:
            return (
                S1_star - 1 + get_mif_len(G2, F2, active_v if active_v in F2 else None)
            )

        else:
            return get_mif_len(
                sg_H, F1, active_v if active_v in F1 else None
            ) + get_mif_len(
                nx.subgraph(G, set(G.nodes) - (H | {v})),
                F - (H | {v}),
                active_v if active_v in F - (H | {v}) else None,
            )

    if len(F) == 0:
        degrees = dict(nx.degree(G))
        v = max(degrees, key=degrees.get)
        new_G = nx.subgraph(G, G.nodes - {v})
        return max(get_mif_len(G, F | {v}, v), get_mif_len(new_G, F, active_v))

    if active_v is None:
        active_v = next(iter(F))

    nb_active = set(nx.neighbors(G, active_v))
    for v in nb_active:
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        d_v = G.degree(v)

        if 5 <= d_v - 1 <= len(gen_nb):
            return max(
                get_mif_len(G, F | {v}, active_v),
                get_mif_len(nx.subgraph(G, set(G.nodes) - {v}), F, active_v),
            )

    for v in nb_active:
        gen_nb = get_generalized_neighbors(G, F, active_v, v)
        if len(gen_nb) == 2:
            if not nx.is_forest(nx.subgraph(G, F | gen_nb)):
                return get_mif_len(G, F | {v}, active_v)
            else:
                return max(
                    get_mif_len(G, F | {v}, active_v),
                    get_mif_len(
                        nx.subgraph(G, set(G.nodes) - {v}),
                        F | gen_nb,
                        active_v,
                    ),
                )

    optimal_v = find_optimal_v(G, F, active_v)
    gen_nb = get_generalized_neighbors(G, F, active_v, optimal_v)
    if len(gen_nb) == 3 and not is_trigger_vertex(G, F, active_v, optimal_v):
        v1, v2, v3 = None, None, None
        # v3 (if possible) not in N(active_v) AND should maximize degree (even if it is in N(active_v) in the end)
        not_in_nb = [x for x in gen_nb if x not in nb_active]
        if len(not_in_nb) > 0:
            v3 = max(not_in_nb, key=lambda x: G.degree(x))
        else:
            v3 = max(gen_nb, key=lambda x: G.degree(x))

        v1, v2 = tuple(gen_nb - {v3})
        if not nx.is_forest(nx.subgraph(G, F | {v1, v2})):
            return max(
                get_mif_len(G, F | {optimal_v}, active_v),
                get_mif_len(
                    nx.subgraph(G, set(G.nodes) - {optimal_v}),
                    F | {v3},
                    active_v,
                ),
            )
        else:
            return max(
                get_mif_len(G, F | {optimal_v}, active_v),
                get_mif_len(
                    nx.subgraph(G, set(G.nodes) - {optimal_v, v3}),
                    F | {v1, v2},
                    active_v,
                ),
                get_mif_len(
                    nx.subgraph(G, set(G.nodes) - {optimal_v}),
                    F | {v3},
                    active_v,
                ),
            )

    else:
        return max(
            get_mif_len(G, F | {optimal_v}, active_v),
            get_mif_len(nx.subgraph(G, set(G.nodes) - {optimal_v}), F, active_v),
        )


def get_mif_len(G, F, active_v):
    if len(F) > 1 and not nx.is_forest(nx.subgraph(G, F)):
        print("Can't reduce cause F is not acyclic")
        return Exception

    new_G = nx.MultiGraph(G)
    # Security measure but shouldn't be necessary, the code shouldn't provoque this case
    new_F = set(F) - set([n for n in F if n not in set(G.nodes)])
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
        for n in set(new_G.nodes) - new_F:
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
            if deg <= 1:
                v = n
                break

        if v is not None:
            S.add(v)
            new_G.remove_node(v)
            new_F.discard(v)
            if active_v == v:
                active_v = None
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

        else:
            break

    if len(new_G.nodes) == 0:
        return len(S)

    else:
        return len(S) + main_procedure(new_G, new_F, active_v)


def get_decycling_number_mif_v3(G):
    if nx.is_forest(G):
        return 0

    return len(G.nodes) - get_mif_len(G, set(), None)


if __name__ == "__main__":
    nt = nx.Graph()
    nt.add_nodes_from(["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"])
    nt.add_edges_from(
        [
            ("V1", "V2"),
            ("V1", "V3"),
            ("V2", "V4"),
            ("V4", "V5"),
            ("V3", "V5"),
            ("V4", "V6"),
            ("V6", "V7"),
            ("V5", "V7"),
            ("V7", "V8"),
            ("V6", "V8"),
        ]
    )

    print(get_decycling_number_mif_v3(nt))
