import networkx as nx
from collections import deque


def find_maximal_2_3_subgraph(og_G):
    G = og_G.copy()
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    nodes_to_visit = deque(G.nodes())

    # Chemin actuel
    stack = []
    in_stack = set()

    # stack[0] est connecté à un noeud de H de degré 2?
    start_connected = False

    while True:
        if len(stack) == 0:
            start_node = None
            while len(nodes_to_visit) > 0:
                n = nodes_to_visit.popleft()
                if G.degree(n) > 0:
                    start_node = n
                    break

            if start_node is None:
                # Plus de noeuds à visiter
                break

            stack = [start_node]
            in_stack = {start_node}
            if H.degree(start_node) == 2:
                start_connected = True
            else:
                start_connected = False

        # Extrémité du chemin actuel
        u = stack[-1]
        parent = stack[-2] if len(stack) > 1 else None
        nb = list(G.neighbors(u))
        valid = []
        edges_to_remove = []
        for n in nb:
            if n != parent:
                if H.degree(n) == 3:
                    edges_to_remove.append((u, n))
                else:
                    valid.append(n)

        if len(edges_to_remove) > 0:
            G.remove_edges_from(edges_to_remove)

        # Si pas de voisins du tout, ou seul voisin est le parent, ou tous les voisins (sans compter le parent)
        # ont un degré 3 dans H -> cul de sac
        if len(valid) == 0:
            if parent is not None:
                if G.has_edge(parent, u):
                    G.remove_edge(parent, u)
            in_stack.remove(stack.pop())

            # Reset si la stack est vide
            if len(stack) == 0:
                start_connected = False
            continue

        v = valid[0]

        # Si v est dans la stack, alors on a un cycle
        if v in in_stack:
            idx = stack.index(v)
            cycle_nodes = stack[idx:]
            edges_to_add = []
            for i in range(len(cycle_nodes) - 1):
                edges_to_add.append((cycle_nodes[i], cycle_nodes[i + 1]))
            edges_to_add.append((u, v))
            H.add_edges_from(edges_to_add)
            G.remove_edges_from(edges_to_add)

            # Si idx == 0, cycle formé par tous les sommets de la stack
            if idx == 0:
                stack = []
                in_stack = set()
                start_connected = False

            else:
                stack = stack[: idx + 1]
                in_stack = set(stack)
                # Le bout de la pile (v) est connecté à H
                # Si le début est connecté, alors c'est un chemin valide à ajouter
                if start_connected:
                    path_edges = []
                    for i in range(len(stack) - 1):
                        path_edges.append((stack[i], stack[i + 1]))
                    H.add_edges_from(path_edges)
                    G.remove_edges_from(path_edges)
                    stack = []
                    in_stack = set()
                    start_connected = False

                else:
                    # La pile est inversée pour essayer d'attacher le début du chemin à H ou un autre cycle
                    stack.reverse()
                    # v était à la fin de la stack, maintenant au début donc le début est connecté à H
                    start_connected = True

            continue

        elif H.degree(v) == 2:
            # v est déjà dans H et à un degré 2 donc la connexion du chemin dans stack est possible à v
            # Si le début est connecté, alors c'est un chemin valide à ajouter
            if start_connected:
                path_edges = []
                for i in range(len(stack) - 1):
                    path_edges.append((stack[i], stack[i + 1]))
                path_edges.append((u, v))
                H.add_edges_from(path_edges)
                G.remove_edges_from(path_edges)
                stack = []
                in_stack = set()
                start_connected = False

            else:
                stack.append(v)
                in_stack.add(v)
                stack.reverse()
                start_connected = True
            continue

        # v est degré 0 ou 1 dans H, on peut l'ajouter au chemin
        else:
            stack.append(v)
            in_stack.add(v)
            continue

    to_remove = [n for n in H.nodes() if H.degree(n) == 0]
    H.remove_nodes_from(to_remove)
    return H


def get_critical_linkpoints(G, H):
    linkpoints = {n for n in H.nodes if H.degree(n) == 2}
    critical_linkpoints = set()
    G_prime = nx.subgraph(G, set(G.nodes) - (set(H.nodes)))
    node_in_component = {}
    for i, comp in enumerate(nx.connected_components(G_prime)):
        for node in comp:
            node_in_component[node] = i

    for n in linkpoints:
        visited_components = set()
        is_critical = False

        for nb in G.neighbors(n):
            if nb in H.nodes:
                continue

            # Vérifier si un cycle comprenant n existe dans G\H avec n étant le sommet du cycle appartenant à H
            if nb in node_in_component:
                comp = node_in_component[nb]
                if comp in visited_components:
                    # Si deux voisins appartiennent à la même composante connexe de G\H, alors il y a un cycle avec n
                    is_critical = True
                    break

                visited_components.add(comp)

        if is_critical:
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
    if nx.is_forest(G):
        return set()

    H = find_maximal_2_3_subgraph(G)
    X = get_critical_linkpoints(G, H)
    Y = {n for n in H.nodes if H.degree(n) >= 3}
    W = get_set_covering_cycles(H, X, Y)
    return X | Y | W


def approx_decycling_number_bar_yehuda(G):
    return len(subG_2_3(G))
