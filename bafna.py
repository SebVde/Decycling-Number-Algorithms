import networkx as nx


def cleanup(G):
    queue = [n for n in G.nodes if G.degree(n) <= 1]

    while len(queue) > 0:
        node = queue.pop()
        if not G.has_node(node):
            continue

        neighbors = list(G.neighbors(node))
        G.remove_node(node)
        for neighbor in neighbors:
            if G.degree(neighbor) <= 1:
                queue.append(neighbor)

    return G


def find_semidisjoint_cycle(G):
    degree_2 = {n for n, d in G.degree() if d == 2}
    if len(degree_2) == 0:
        return None

    sg_d2 = nx.subgraph(G, degree_2)
    visited = set()

    for node in degree_2:
        if node in visited:
            continue

        component = list(nx.node_connected_component(sg_d2, node))
        visited.update(component)
        sg_comp = nx.subgraph(G, component)

        # If in the component the number of edges equals the number of nodes, then it's a cycle
        if 0 < len(component) == sg_comp.number_of_edges():
            return component

        # Else, we need to check if it's a path with endpoints connected to a common junction node
        endpoints = [n for n in component if sg_comp.degree(n) <= 1]
        if len(endpoints) == 2:
            ep1, ep2 = endpoints
            nb_1 = {n for n in G.neighbors(ep1) if n not in component}
            nb_2 = {n for n in G.neighbors(ep2) if n not in component}
            common_junctions = {n for n in nb_1 & nb_2 if G.degree(n) > 2}

            if len(common_junctions) > 0:
                return component + [list(common_junctions)[0]]

    return None


def get_fvs(og_G):
    F = set()
    stack = []
    G = og_G.copy()
    weights = {n: 1.0 for n in G.nodes()}

    while len(G.nodes) > 0:
        sd_cycle = find_semidisjoint_cycle(G)
        if sd_cycle is not None:
            for node in sd_cycle:
                weights[node] -= min(weights[node] for node in sd_cycle)

        else:
            gamma = min(weights[node] / (G.degree(node) - 1) for node in G.nodes)
            for node in G.nodes:
                weights[node] -= gamma * (G.degree(node) - 1)

        to_remove = [node for node in G.nodes if weights[node] <= 1e-9]
        F.update(to_remove)
        G.remove_nodes_from(to_remove)
        stack.extend(to_remove)
        G = cleanup(G)

    while len(stack) > 0:
        node = stack.pop()
        sg = nx.subgraph_view(og_G, filter_node=lambda n: n not in F or n == node)
        # With this check, it should ensure that each call for nx.is_forest is in O(|V|) time since |E| < |V|
        if sg.number_of_edges() <= sg.number_of_nodes() - 1 and nx.is_forest(sg):
            F.remove(node)

    return F


def approx_decycling_number_bafna(G):
    clean_G = cleanup(G.copy())
    fvs = get_fvs(clean_G)
    return len(fvs)
