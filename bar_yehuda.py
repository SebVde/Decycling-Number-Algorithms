import networkx as nx
from collections import deque


def find_maximal_2_3_subgraph(og_G):
    """
    Finds a maximal subgraph of the input graph where all nodes have degree 2 or 3.
    This function iteratively constructs a subgraph by identifying cycles and paths
    that can be added while maintaining the degree constraints.

    Args:
        og_G (nx.Graph): The input graph.

    Returns:
        nx.Graph: A maximal subgraph of the input graph where all nodes have degree 2 or 3.
    """
    G = og_G.copy()
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    nodes_to_visit = deque(G.nodes())

    # Current path
    stack = []
    in_stack = set()

    # Is stack[0] connected to H by a degree 2 node?
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
                # No more nodes to visit
                break

            stack = [start_node]
            in_stack = {start_node}
            if H.degree(start_node) == 2:
                start_connected = True
            else:
                start_connected = False

        # Extract the extremity of the path
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

        # If no valid neighbors, backtrack
        if len(valid) == 0:
            if parent is not None:
                if G.has_edge(parent, u):
                    G.remove_edge(parent, u)
            in_stack.remove(stack.pop())

            # Reset if stack is empty
            if len(stack) == 0:
                start_connected = False
            continue

        v = valid[0]

        # If v is already in stack -> cycle detected
        if v in in_stack:
            idx = stack.index(v)
            cycle_nodes = stack[idx:]
            edges_to_add = []
            for i in range(len(cycle_nodes) - 1):
                edges_to_add.append((cycle_nodes[i], cycle_nodes[i + 1]))
            edges_to_add.append((u, v))
            H.add_edges_from(edges_to_add)
            G.remove_edges_from(edges_to_add)

            # Handle cycle cases
            if idx == 0:
                stack = []
                in_stack = set()
                start_connected = False
            else:
                stack = stack[: idx + 1]
                in_stack = set(stack)
                # The end of the stack (v) is connected to H
                # If the start is connected, then it's a valid path to add
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
                    # The stack is reversed to try to connect the start of the path to H or another cycle
                    stack.reverse()
                    # v was at the end of the stack, now at the beginning so the start of the stack is connected to H
                    start_connected = True

            continue

        elif H.degree(v) == 2:
            # Connect path to v if start is connected
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

        # v has degree 0 or 1 in H, we can add it to the path
        else:
            stack.append(v)
            in_stack.add(v)
            continue

    to_remove = [n for n in H.nodes() if H.degree(n) == 0]
    H.remove_nodes_from(to_remove)
    return H


def get_critical_linkpoints(G, H):
    """
    Identifies critical linkpoints in the subgraph.
    A critical linkpoint v is a node of degree 2 in H such that removing
    all the nodes of H except v from G creates a cycle including v.

    Args:
        G (nx.Graph): The original graph.
        H (nx.Graph): The maximal subgraph with degree 2 or 3 nodes.

    Returns:
        set: The set of critical linkpoints.
    """
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

            # Verify if a cycle including n exists in G\H with n being the vertex of the cycle belonging to H
            if nb in node_in_component:
                comp = node_in_component[nb]
                if comp in visited_components:
                    # If two neighbors belong to the same connected component of G\H, then there is a cycle with n
                    is_critical = True
                    break

                visited_components.add(comp)

        if is_critical:
            critical_linkpoints.add(n)

    return critical_linkpoints


def is_cycle(G):
    """
    Checks if the graph is a cycle.

    Args:
        G (nx.Graph): The input graph.

    Returns:
        bool: True if the graph is a cycle, False otherwise.
    """
    return all(G.degree(n) == 2 for n in G.nodes)


def get_set_covering_cycles(H, X, Y):
    """
    Finds a set of nodes that cover all cycles in the subgraph excluding critical linkpoints and high-degree nodes.

    Args:
        H (nx.Graph): The subgraph with degree 2 or 3 nodes.
        X: Set of critical linkpoints.
        Y: Nodes with degree = 3 in H.

    Returns:
        set: A set of nodes covering all cycles in the subgraph.
    """
    sg = nx.subgraph(H, set(H.nodes) - X - Y)
    cover_set = set()

    for comp in nx.connected_components(sg):
        comp_sg = nx.subgraph(sg, comp)
        if is_cycle(comp_sg):
            cover_set.add(next(iter(comp)))

    return cover_set


def subG_2_3(G):
    """
    Approximates a minimal FVS for G using Bar-Yehuda's algorithm.

    Args:
        G (nx.Graph): The input graph.

    Returns:
        set: An approximated minimal FVS for G.
    """
    if nx.is_forest(G):
        return set()

    H = find_maximal_2_3_subgraph(G)
    X = get_critical_linkpoints(G, H)
    Y = {n for n in H.nodes if H.degree(n) >= 3}
    W = get_set_covering_cycles(H, X, Y)
    return X | Y | W


def approx_decycling_number_bar_yehuda(G):
    """
    Approximates the decycling number of the graph using Bar-Yehuda's algorithm.

    Args:
        G (nx.Graph): The input graph.

    Returns:
        int: The approximated decycling number.
    """
    return len(subG_2_3(G))
