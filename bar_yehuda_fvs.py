import networkx as nx
from collections import deque


def dfs_construct(G, H, v, visited, edges_visited):
    visited.add(v)
    for u in G.neighbors(v):
        edge = (min(v, u), max(v, u))
        if edge in edges_visited:
            continue

        edges_visited.add(edge)
        deg_u = H.degree[u] if u in H else 0
        deg_v = H.degree[v] if v in H else 0

        if deg_u < 3 and deg_v < 3:
            H.add_edge(u, v)
        if u not in visited:
            dfs_construct(G, H, u, visited, edges_visited)


def find_maximal_2_3_subgraph(G):
    H = nx.Graph()
    visited = set()
    edges_visited = set()
    for node in G.nodes:
        if node not in visited:
            dfs_construct(G, H, node, visited, edges_visited)

    to_remove = {n for n in H.nodes if H.degree(n) < 2}
    while len(to_remove) > 0:
        v = to_remove.pop()
        if v not in H:
            continue

        nb = set(H.neighbors(v))
        H.remove_node(v)
        for u in nb:
            if u in H and H.degree[u] < 2:
                to_remove.add(u)

    return H


def cycle_exists_with_node(G, n):
    queue = deque()
    queue.append((n, None))
    visited = set()
    visited.add(n)

    while queue:
        current, parent = queue.popleft()
        for nb in G.neighbors(current):
            if nb == parent:
                continue

            if nb == n:
                return True

            if nb not in visited:
                visited.add(nb)
                queue.append((nb, current))

    return False


def get_critical_linkpoints(G, H):
    linkpoints = {n for n in H.nodes if H.degree(n) == 2}
    critical_linkpoints = set()

    for n in linkpoints:
        sg = nx.subgraph(G, set(G.nodes) - (set(H.nodes) - {n}))
        if cycle_exists_with_node(sg, n):
            critical_linkpoints.add(n)

    return critical_linkpoints


def is_cycle(G):
    return all(G.degree(n) == 2 for n in G.nodes)


def get_set_covering_cycles(H, X, Y):
    sg = nx.subgraph(H, set(H.nodes) - X - Y)
    cover_set = set()

    for comp in nx.connected_components(sg):
        comp_sg = nx.subgraph(sg, comp)
        if is_cycle(comp_sg):
            cover_set.add(next(iter(comp)))

    return cover_set


def subG_2_3(G):
    # Complexité pas linéaire à cause de get_critical_linkpoints
    if nx.is_forest(G):
        return set()

    H = find_maximal_2_3_subgraph(G)
    X = get_critical_linkpoints(G, H)
    Y = {n for n in H.nodes if H.degree(n) >= 3}
    W = get_set_covering_cycles(H, X, Y)
    return X | Y | W


def get_decycling_number_bar_yehuda(G):
    return len(subG_2_3(G))
