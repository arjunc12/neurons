import numpy as np
import networkx as nx

from sys import argv

def makes_cycle(u, v, node_to_forest):
    f1 = node_to_forest[u]
    f2 = node_to_forest[v]
    return f1 == f2
 
def combine_forests(u, v, node_to_forest, forest_to_nodes, forest_to_edges):
    f1 = node_to_forest[u]
    f2 = node_to_forest[v]
    assert f1 != f2
    new_forest = min(f1, f2)
    old_forest = max(f1, f2)
    update_nodes = forest_to_nodes[old_forest]
    update_edges = forest_to_edges[old_forest]
    for node in update_nodes:
        node_to_forest[node] = new_forest
    forest_to_nodes[new_forest] += update_nodes
    forest_to_edges[new_forest] += update_edges
    forest_to_edges[new_forest].append((u, v))
    del forest_to_nodes[old_forest]
    del forest_to_edges[old_forest]
    
def kruskal(nodes, edges):
    node_to_forest = {}
    forest_to_nodes = {}
    forest_to_edges = {}
    
    for i, u in enumerate(nodes):
        node_to_forest[u] = i
        forest_to_nodes[i] = [u]
        forest_to_edges[i] = []
        
    for u, v in edges:
        if not makes_cycle(u, v, node_to_forest):
            combine_forests(u, v, node_to_forest, forest_to_nodes, forest_to_edges)
            
    return sum(forest_to_edges.values(), [])

def mst_cost(G):
    total_length = 0
    for u, v in G.edges_iter():
        total_length += G[u][v]['length']
    return total_length

def satellite_cost(G, root):
    total_cost = 0
    has_root_path = False
    for u in G.nodes_iter():
        if u == root:
            continue
        if nx.has_path(G, root, u):
            total_cost += nx.shortest_path_length(G, root, u, weight='length')
            has_root_path = True
    if has_root_path:
        return total_cost
    else:
        return float("inf")

def pareto_kruskal(G, root, alpha):
    pareto_mst = nx.Graph()
    node_to_forest = {}
    forest_to_nodes = {}
    forest_to_edges = {}

    for i, u in enumerate(G.nodes()):
        node_to_forest[u] = i
        forest_to_nodes[i] = [u]
        forest_to_edges[i] = []
    
    pareto_mst.add_node(root)

    added_edges = set()
    ignore_edges = set()
    candidate_edges = G.edges()
    while len(forest_to_nodes) > 1:
        best_edge = None
        best_cost = float('inf')
        for u, v in candidate_edges:
            if (u, v) in added_edges or (u, v) in ignore_edges:
                continue
            if makes_cycle(u, v, node_to_forest):
                ignore_edges.add((u, v))
                continue
            pareto_mst.add_edge(u, v)
            pareto_mst[u][v]['length'] = G[u][v]['length']
            mcost = alpha * mst_cost(pareto_mst)

            scost = (1 - alpha) * satellite_cost(pareto_mst, root)
            cost = mcost + scost
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
            pareto_mst.remove_edge(u, v)
        if best_edge == None:
            break
        u, v = best_edge
        added_edges.add((u, v))
        pareto_mst.add_edge(u, v)
        pareto_mst[u][v]['length'] = G[u][v]['length']
        combine_forests(u, v, node_to_forest, forest_to_nodes, forest_to_edges)
    return pareto_mst

def point_dist(p1, p2):
    assert len(p1) == len(p2)
    sq_dist = 0
    for i in xrange(len(p1)):
        x1, x2 = p1[i], p2[i]
        sq_dist += (x1 - x2) ** 2
    return sq_dist ** 0.5

def pareto_mst(points, root, alpha):
    G = nx.Graph()
    for i in xrange(len(points)):
        for j in xrange(len(points)):
            if i == j:
                continue
            p1, p2 = points[i], points[j]
            G.add_edge(p1, p2)
            length = point_dist(p1, p2)
            G[p1][p2]['length'] = length
    mst = pareto_kruskal(G, root, alpha)
    return mst

if __name__ == '__main__':
    points = [(0, 0), (1, 1), (1, 1.1), (0, 0.1), (2, 2)]
    root = (0, 0)
    alpha = float(argv[1])
    mst = pareto_mst(points, root, alpha)
    print mst.nodes()
    print mst.edges()
