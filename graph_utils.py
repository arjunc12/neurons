import networkx as nx
from itertools import combinations
from dist_functions import point_dist

def complete_graph(G):
    H = G.copy()
    H.remove_edges_from(H.edges())
    for u, v in combinations(H.nodes(), 2):
        H.add_edge(u, v)
        H[u][v]['length'] = point_dist(H.node[u]['coord'], H.node[v]['coord'])
    return H

def is_tree(G):
    if G.number_of_edges() != G.number_of_nodes() - 1:
        print "wrong number of edges"
        print "edges", G.number_of_edges(), "nodes", G.number_of_nodes()
        return False
    if not nx.is_connected(G):
        print "not connected"
        return False
    return True

def check_dists(G):
    assert 'root' in G.graph
    root = G.graph['root']
    assert 'coord' in G.node[root]
    root_coord = G.node[root]['coord']
    for u in G.nodes():
        assert 'droot' in G.node[u]
        droot = G.node[u]['droot']
        assert 'coord' in G.node[u]
        coord = G.node[u]['coord']
        min_dist = point_dist(coord, root_coord)
        if droot < min_dist:
            print "bad node", u
            print "droot", droot, "min_dist", min_dist
            print type(droot), type(min_dist)
            print droot - min_dist
        assert droot >= min_dist

def dist_error(G):
    assert 'root' in G.graph
    root = G.graph['root']
    assert 'coord' in G.node[root]
    root_coord = G.node[root]['coord']
    error = 0
    for u in G.nodes():
        assert 'droot' in G.node[u]
        droot = G.node[u]['droot']
        assert 'coord' in G.node[u]
        coord = G.node[u]['coord']
        min_dist = point_dist(coord, root_coord)
        if droot < min_dist:
            error += min_dist - droot

    return error

def scost_error(G):
    scost = 0
    scost_error = 0
    
    droot = {}
    assert 'root' in G.graph
    root = G.graph['root']
    droot[root] = 0
    assert 'coord' in G.node[root]
    root_coord = G.node[root]['coord']
    
    parent = {}
    parent[root] = None

    queue = [root]
    visited = set()
    while len(queue) > 0:
        curr = queue.pop(0)
        
        if curr in visited:
            print "not a tree; graph has cycles"
            assert False
        
        visited.add(curr)

        for child in G.neighbors(curr):
            if child != parent[curr]:
                assert 'length' in G[curr][child]
                length = G[curr][child]['length']
                
                child_droot = length + droot[curr]
                scost += child_droot
                assert 'coord' in G.node[child]
                coord = G.node[child]['coord']
                min_dist = point_dist(root_coord, coord)
                if child_droot < min_dist:
                    error = min_dist - child_droot
                    scost_error += error
                
                droot[child] = child_droot
                parent[child] = curr
                queue.append(child)

    if len(visited) < G.number_of_nodes():
        scost = float("inf") # graph is not connected

    return scost, scost_error

def points_to_graph(points):
    point_graph = nx.Graph()
    for p1, p2 in combinations(points, 2):
        point_graph.add_edge(p1, p2)
        length = point_dist(p1, p2)
        point_graph[p1][p2]['length'] = length
    return point_graph

def non_continue_subgraph(G):
    H = G.copy()
    root = H.graph['root']
    for u in G.nodes_iter():
        if G.node[u]['label'] == 'continue':
            H.remove_node(u)
    H = complete_graph(H)
    return H

def graph_to_points(G, root, P2Coord, cont=False):
    points = []
    for u in G.nodes_iter():
        degree = len(G.neighbors(u))
        if (u == root) or (u != root and degree != 2) or cont:
            points.append(P2Coord[u])
    return points

def init_lengths(G):
    for u, v in G.edges_iter():
        p1, p2 = G.node[u]['coord'], G.node[v]['coord']
        G[u][v]['length'] = point_dist(p1, p2)
