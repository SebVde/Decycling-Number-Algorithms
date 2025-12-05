import networkx as nx


def dfs_remove_below_2(H, v):
    # Supprime récursivement les sommets de degré inférieur à 2
    # Pas de visited car vu qu'on supprime des noeuds, certains noeuds peuvent changer de degré
    if (v in H) and (H.degree(v) < 2):
        nb = set(H.neighbors(v))
        H.remove_node(v)

        for u in nb:
            dfs_remove_below_2(H, u)


def find_maximal_2_3_subgraph(G):
    # Sous-graphe maximal de G avec que des degrés entre 2 et 3 dans ce sous-graphe
    H = G.copy()

    while any(H.degree(node) > 3 for node in H.nodes):
        node = next(n for n in H.nodes if H.degree(n) > 3)
        node_deg = H.degree(node)
        if node_deg > 3:
            m = node_deg - 3
            nb = sorted(
                list(H.neighbors(node)), key=lambda x: H.degree(x), reverse=True
            )
            for i in range(m):
                H.remove_edge(node, nb[i])

    for node in set(H.nodes):
        dfs_remove_below_2(H, node)
    return H


def subG_2_3(G):
    if nx.is_forest(G):
        return set()

    H = find_maximal_2_3_subgraph(G)
    pass


def get_decycling_number_bar_yehuda(G):
    return len(subG_2_3(G))


if __name__ == "__main__":
    G = nx.Graph()
    edges = [
        (1, 2),
        (2, 3),
        (3, 4),
        (2, 4),
        (3, 5),
        (5, 6),
        (6, 7),
        (7, 5),
    ]
    G.add_edges_from(edges)

    # ajouter des sommets de degré 1 (feuilles) au graphe G existant
    G.add_edge(1, 8)  # 8 devient un sommet de degré 1
    G.add_edge(6, 9)  # 9 devient un sommet de degré 1
    G.add_edge(7, 10)  # 10 devient un sommet de degré 1

    print("Graph G:")
    print("Nodes:", G.nodes())
    print("Edges:", G.edges())
    #
    # decycling_number = get_decycling_number_bar_yehuda(G)
    # print("Decycling number (Bar-Yehuda method):", decycling_number)
    H = find_maximal_2_3_subgraph(G)
    print(H)
