import networkx as nx
import random


class CycleDetected(Exception):
    def __init__(self, message):
        super().__init__(message)


# Algorithm taken from the article "Exact Computation of Maximum Induced Forest"


def get_cnf(G, T):
    """
Computes the set of conflicting vertices in G that are not in T but have at least 2 neighbors in a same connected
component of the subgraph of G induced by T.
    :param G: a networkx graph
    :param T: a subset of vertices in G
    :return: a set of vertices representing the conflicting vertices
    """

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
                break  # Because v only needs to have 2 neighbors in one of the components to be considered conflicting

    return cnf


def get_bnd(G, T):
    """
Computes the boundary vertices of T in G ie the vertices that are not in T but have exactly one neighbor in T.
    :param G: a networkx graph
    :param T: a subset of vertices in G
    :return: a set of vertices representing the boundary vertices of T in G
    """

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
    """
Computes the subset of vertices in K that are K-connected to v.
To be K-connected to v, a vertex u in K must have a path to v in G.
    :param G: a networkx graph
    :param K: a subset of vertices in G
    :param v: a vertex in G
    :return: a set of vertices representing the K-connected vertices to v
    """

    result = set()
    for u in K:
        if nx.has_path(G, u, v):
            result.add(u)

    return result


def T_update(G, T, K, v):
    S = get_K_connected(G, K, v)
    new_T = T | S | {v}  # Update T with the K-connected vertices and v
    new_K = K - S - {v}  # Update K by removing the K-connected vertices and v
    sg_nodes = G.nodes - get_cnf(G, new_T)
    sg = nx.subgraph(G, sg_nodes)  # Create a subgraph induced by the nodes that are not conflicting with the new T
    # All these conflicting vertices are removed because they would create cycles in the subgraph
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
        return None  # If the subgraph contains cycles


def find_mif_len(G, T, K):
    """
    Searches for the length of a (T union K)-MIF of the graph G.
    :param G: a networkx graph
    :param T: a subset of vertices in G
    :param K: a subset of vertices in G
    :return: the length of a (T union K)-MIF of G
    """

    if T | K == set(G.nodes):
        return len(T | K)

    else:
        v = None
        bnd = get_bnd(G, T)
        if len(bnd) != 0:
            v = random.choice(list(bnd))
        else:
            v = random.choice(list(G.nodes - (T | K)))

        new_G, new_T, new_K = T_update(G, T, K, v)
        new_S = find_mif_len(new_G, new_T, new_K)

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
                result = K_update(sg, T, K, W)
                if result is not None:
                    return max(new_S, find_mif_len(*result))
                else:
                    return new_S
            case _:
                return max(
                    new_S,
                    find_mif_len(nx.subgraph(G, G.nodes - {v}), T, K),
                )


def get_mif_len(G, K):
    """
Main function to get the length of a K-MIF of a graph G.
Given a subset K of V(G), a K-MIF is a largest superset of K such that the subgraph of G induced by K is acyclic.
    :param G: a networkx graph
    :param K: a subset of vertices in G
    :raises CycleDetected: if the subgraph induced by K contains cycles
    :return: the length of a K-MIF of G
    """

    if len(K) == 0 or nx.is_forest(nx.subgraph(G, K)):
        sg_nodes = G.nodes - get_cnf(G, K)
        sg = nx.subgraph(G, sg_nodes)
        result = find_mif_len(sg, set(), K)
        return result

    else:
        raise CycleDetected(
            "No K-MIF of G can be found because the subgraph induced by K contains cycles."
        )


def get_decycling_number_mif(G):
    return len(G.nodes) - get_mif_len(G, set())
