import networkx as nx

def normalize_cost(cost, opt_cost):
    return 1 - (opt_cost / cost)

def graph_costs(G):
    scost = 0
    mcost = 0
    
    droot = {}
    root = G.graph['root']
    droot[root] = 0
    
    parent = {}
    parent[root] = None

    queue = [root]
    visited = set()
    while len(queue) > 0:
        curr = queue.pop(0)
        
        if curr in visited:
            return float("inf"), float("inf")
        
        visited.add(curr)
        for child in G.neighbors(curr):
            if child != parent[curr]:
                length = G[curr][child]['length']
                mcost += length
                child_droot = length + droot[curr]
                droot[child] = child_droot
                scost += child_droot
                parent[child] = curr
                queue.append(child)

    if len(visited) < G.number_of_nodes():
        scost = float("inf") # graph is not connected

    return mcost, scost

def best_mst_cost(G):
    sorted_edges = sorted(G.edges(), key= lambda x : G[x[0]][x[1]]['length'])
    mst_edges = kruskal(G.nodes(), sorted_edges)
    best_cost = 0
    for u, v in mst_edges:
        best_cost += G[u][v]['length']
    return best_cost

def mst_cost(G):
    total_length = 0
    for u, v in G.edges_iter():
        total_length += G[u][v]['length']
    return total_length

def best_satellite_cost(G):
    best_cost = 0
    root = G.graph['root']
    root_coord = G.node[root]['coord']
    for u in G.nodes_iter():
        best_cost += point_dist(root_coord, G.node[u]['coord'])
    return best_cost

def satellite_cost(G):
    root = G.graph['root']
    total_cost = 0
    for u in G.nodes_iter():
        if u == root:
            continue
        if nx.has_path(G, root, u):
            total_cost += nx.shortest_path_length(G, root, u, weight='length')
        else:
            total_cost += float("inf")
    return total_cost
