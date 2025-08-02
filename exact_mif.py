import networkx as nx
import random


class CycleDetected(Exception):
    def __init__(self, message):
        super().__init__(message)


# Algorithm taken from the article "Exact Computation of Maximum Induced Forest"


def get_cnf(G, T):
    cnf = set()
    components = list(nx.connected_components(nx.subgraph(G, T)))

    for v in G.nodes():
        if v in T:
            continue

        # For each connected component from T, check if v has 2 or more neighbors in the component
        for comp in components:
            neighbors_in_comp = set(nx.neighbors(G, v)) & comp  # & = intersection
            if len(neighbors_in_comp) >= 2:
                cnf.add(v)
                break

    return cnf


def get_bnd(G, T):
    bnd = set()

    for v in G.nodes():
        if v in T:
            continue

        neighbors_to_T = set(nx.neighbors(G, v)) & T
        if len(neighbors_to_T) == 1:
            bnd.add(v)
            break

    return bnd


def get_K_connected(G, K, v):
    # Computes the subset of vertices in K that are K-connected to v
    result = set()
    for u in K:
        if nx.has_path(G, u, v):
            result.add(u)

    return result


def T_update(G, T, K, v):
    S = get_K_connected(G, K, v)
    new_T = T | S | {v}
    new_K = K - S - {v}
    sg_nodes = G.nodes - get_cnf(G, new_T)
    sg = nx.subgraph(G, sg_nodes)
    return sg, new_T, new_K


def K_update(G, T, K, W):
    sg = nx.subgraph(G, T | K | W)

    if nx.is_forest(sg):
        bnd = get_bnd(G, T)
        if len(W & bnd) == 0:
            sg_nodes = G.nodes - get_cnf(G, K | W)
            return nx.subgraph(G, sg_nodes), T, K | W

        else:
            v1 = random.choice(list(W & bnd))
            new_G, new_T, new_K = T_update(G, T, K, v1)
            v2 = W - {v1}

            if v2 in get_bnd(new_G, new_T):
                return T_update(new_G, new_T, new_K, v2)
            else:
                sg_nodes = new_G.nodes - get_cnf(new_G, new_K | {v2})
                return nx.subgraph(new_G, sg_nodes), new_T, new_K | {v2}

    else:
        return None  # If the subgraph is not a forest, return Fail


def find_mif(G, T, K):
    if T | K == set(G.nodes):
        return T | K

    else:
        v = None
        bnd = get_bnd(G, T)
        if len(bnd) != 0:
            v = random.choice(list(bnd))
        else:
            v = random.choice(list(G.nodes - (T | K)))

        new_G, new_T, new_K = T_update(G, T, K, v)
        new_S = find_mif(new_G, new_T, new_K)

        W = get_K_connected(
            nx.subgraph_view(
                G, filter_node=lambda n: n in G.nodes and n not in T | K | {v}
            ),
            K,
            v,
        )

        match len(W):
            case 0, 1:
                return new_S
            case 2:
                sg = nx.subgraph(G, G.nodes - {v})
                result = K_update(sg, new_T, new_K, W)
                if result is not None:
                    return max(new_S, find_mif(*result), key=lambda x: len(x))
                else:
                    return new_S
            case _:
                return max(
                    new_S,
                    find_mif(nx.subgraph(G, G.nodes - {v}), T, K),
                    key=lambda x: len(x),
                )


def main_mif(G, K):
    if len(K) == 0 or nx.is_forest(nx.subgraph(G, K)):
        sg_nodes = G.nodes - get_cnf(G, K)
        sg = nx.subgraph(G, sg_nodes)
        result = sorted(find_mif(sg, set(), K))

        print("MIF found:", result)
        return result

    else:
        raise CycleDetected(
            "The subgraph induced by the given K should NOT contain cycles"
        )
