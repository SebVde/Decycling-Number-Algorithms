import networkx as nx
import heapq
import itertools
from collections import deque


class PriorityQueue:
    def __init__(self):
        self.pq = []
        self.entry_finder = {}
        self.counter = itertools.count()
        self.REMOVED = "<removed-task>"

    def insert(self, item, priority):
        if item in self.entry_finder:
            self.delete(item)
        count = next(self.counter)
        entry = [priority, count, item]
        self.entry_finder[item] = entry
        heapq.heappush(self.pq, entry)

    def delete(self, item):
        entry = self.entry_finder.pop(item, None)
        if entry:
            entry[-1] = self.REMOVED

    def extract_min(self):
        while self.pq:
            priority, count, item = heapq.heappop(self.pq)
            if item is not self.REMOVED:
                del self.entry_finder[item]
                return item, priority
        raise KeyError("pop from an empty priority queue")

    def __len__(self):
        return len(self.entry_finder)


def cleanup(G):
    while True:
        to_remove = [node for node in G.nodes if G.degree(node) <= 1]
        if len(to_remove) == 0:
            break
        G.remove_nodes_from(to_remove)

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


def get_current_weight(G, n, node_data, total_gamma_clean):
    base_w = node_data[n]["weight"]
    last_g = node_data[n]["last_seen_gamma"]
    deg = G.degree(n)
    # w_actuel = w_base - (gamma accumulé - gamma au dernier update) * (degré actuel - 1)
    return base_w - (total_gamma_clean - last_g) * (deg - 1)


# Fonction pour "fixer" le poids actuel (commit) avant un changement de degré
def commit_weight(G, n, node_data, total_gamma_clean):
    curr_w = get_current_weight(G, n, node_data, total_gamma_clean)
    node_data[n]["weight"] = curr_w
    node_data[n]["last_seen_gamma"] = total_gamma_clean
    return curr_w


def get_fvs(og_G):
    F = set()
    G = og_G.copy()
    stack = []

    # pour accumuler les gamma soustraits lorsque le graphe n'a pas de cycle semi-disjoint
    total_gamma_clean = 0.0

    # stocker la valeur qu'avait total_gamma_clean la dernière fois qu'on a touché au noeud n
    node_data = {}
    for n in G.nodes:
        node_data[n] = {"weight": 1.0, "last_seen_gamma": 0.0}

    pq = PriorityQueue()
    for n in G.nodes:
        if G.degree(n) > 1:
            # clé = w(u)/(d(u)-1) + total_gamma_clean (0)
            ratio = node_data[n]["weight"] / (G.degree(n) - 1)
            pq.insert(n, ratio)

    while len(G.nodes) > 0:
        to_remove = []
        sd_cycle = find_semidisjoint_cycle(G)
        if sd_cycle is not None:
            min_w = 2
            for n in sd_cycle:
                w = get_current_weight(G, n, node_data, total_gamma_clean)
                if w < min_w:
                    min_w = w

            gamma = min_w
            for n in sd_cycle:
                commit_weight(G, n, node_data, total_gamma_clean)
                node_data[n]["weight"] -= gamma

                if node_data[n]["weight"] <= 1e-9:
                    to_remove.append(n)
                else:
                    new_ratio = node_data[n]["weight"] / (G.degree(n) - 1)
                    pq.insert(n, new_ratio + total_gamma_clean)

        else:
            if len(pq) == 0:
                break

            try:
                u, ratio = pq.extract_min()
            except KeyError:
                break

            current_ratio = ratio - total_gamma_clean
            gamma = current_ratio
            total_gamma_clean += gamma

            node_data[u][
                "weight"
            ] = 0  # car pour le sommet à min ratio, w(u) = w(u) - (w(u)/(d(u)-1)) * (d(u)-1) = w(u)-w(u)
            node_data[u]["last_seen_gamma"] = total_gamma_clean
            to_remove.append(u)

        queue = deque(to_remove)
        while queue:
            node = queue.popleft()
            if not G.has_node(node):
                continue

            nb = set(G.neighbors(node))
            F.add(node)
            stack.append(node)
            pq.delete(node)
            G.remove_node(node)

            for u in nb:
                u_curr_weight = node_data[u]["weight"] - (
                    total_gamma_clean - node_data[u]["last_seen_gamma"]
                ) * G.degree(
                    u
                )  # degré+1 (anncien degré) - 1

                node_data[u]["weight"] = u_curr_weight
                node_data[u]["last_seen_gamma"] = total_gamma_clean

                if G.degree(u) <= 1:
                    if u not in queue:
                        queue.append(u)
                else:
                    # ratio = poids actuel / (degré actuel - 1)
                    new_ratio = u_curr_weight / (G.degree(u) - 1)
                    pq.insert(u, new_ratio + total_gamma_clean)

    while len(stack) > 0:
        node = stack.pop()
        sg = nx.subgraph_view(og_G, filter_node=lambda n: n not in F or n == node)
        # With this check, it should ensure that each call for nx.is_forest is in O(|V|) time since |E| < |V|
        if sg.number_of_edges() <= sg.number_of_nodes() - 1 and nx.is_forest(sg):
            F.remove(node)

    return F


def get_decycling_number_2_approx(G):
    clean_G = cleanup(G.copy())
    fvs = get_fvs(clean_G)
    return len(fvs)
