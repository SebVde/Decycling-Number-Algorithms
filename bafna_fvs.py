import networkx as nx


def cleanup(G):
    new_G = G.copy()
    while True:
        to_remove = [node for node in new_G.nodes if new_G.degree(node) <= 1]
        if len(to_remove) == 0:
            break
        new_G.remove_nodes_from(to_remove)

    return new_G


def find_semidisjoint_cycle(G):
    for cycle in nx.cycle_basis(G):
        cnt = 0
        for node in cycle:
            if G.degree(node) != 2:
                cnt += 1
            if cnt > 1:
                break
        if cnt <= 1:
            return cycle
    return None


def get_fvs(og_G):
    F = set()
    i = 0
    stack = []
    G = og_G.copy()

    while len(G.nodes) > 0:
        i += 1
        sd_cycle = find_semidisjoint_cycle(G)
        if sd_cycle is not None:
            for node in sd_cycle:
                G.node[node]["weight"] = 0.0

        else:
            gamma = min(1 / (G.degree(node) - 1) for node in G.nodes)
            for node in G.nodes:
                G.nodes[node]["weight"] -= gamma * (G.degree(node) - 1)

        to_remove = [node for node in G.nodes if G.nodes[node]["weight"] == 0.0]
        F.update(to_remove)
        G.remove_nodes_from(to_remove)
        stack.extend(to_remove)
        G = cleanup(G)

    while len(stack) > 0:
        node = stack.pop()
        sg = nx.subgraph(og_G, og_G.nodes - (F - {node}))
        if nx.is_forest(sg):
            F.remove(node)

    return F


def get_decycling_number_2_approx(G):
    clean_G = cleanup(G)
    for node in clean_G.nodes:
        clean_G.nodes[node]["weight"] = 1.0

    fvs = get_fvs(clean_G)
    return len(fvs)
