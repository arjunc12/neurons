import numpy as np
import networkx as nx
from sys import argv
import matplotlib as mpl
mpl.use('agg')
import pylab

def get_neuron_graph(filename):
    #for arbor_type in ["2","3","4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
    for arbor_type in ["3"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
        
        # Reads in 3D arborization.
        G = nx.Graph()
        P2Coord = {} # u-> (x,y)
        root = -1
        with open(filename) as f:
            for line in f:
                if line.startswith("#"): continue

                cols = line.strip().split()
                assert len(cols) == 7

                if not (cols[1] == "1" or cols[1] == arbor_type): continue
                
                if cols[6] == "-1":
                    if root != -1: assert False # duplicate root.

                    root = int(cols[0])
                                        
                    assert root not in P2Coord and root not in G
                    G.add_node(root)
                    
                    if dim == "3D":                
                        P2Coord[root] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        P2Coord[root] = (float(cols[2]),float(cols[3]))
                    else:
                        assert False


                else:
                    u,v = int(cols[6]),int(cols[0])
                    assert u in G and u in P2Coord
                    assert not G.has_edge(u,v)
            
                    if dim == "3D":
                        P2Coord[v] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        P2Coord[v] = (float(cols[2]),float(cols[3]))
                    else: 
                        assert False

                    G.add_edge(u,v)

                    #if UTIL.euclidean_dist(u,v,P2Coord) > 10:
                    #    print u,v
        
        
    

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
    for u in G.nodes_iter():
        if u == root:
            continue
        if nx.has_path(G, root, u):
            total_cost += nx.shortest_path_length(G, root, u, weight='length')
        else:
            total_cost += float("inf")
    return total_cost

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
            mst_copy = pareto_mst.copy()
            mst_copy.add_edge(u, v)
            mst_copy[u][v]['length'] = G[u][v]['length']
            mcost = alpha * mst_cost(mst_copy)

            scost = (1 - alpha) * satellite_cost(mst_copy, root)
            cost = mcost + scost
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
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

def pareto_plot(points, root):
    mcosts = []
    scosts = []
    for alpha in np.arange(0, 1.05, 0.05):
        mst = pareto_mst(points, root, alpha)
        mcost = mst_cost(mst)
        scost = satellite_cost(mst, root)
        
        mcosts.append(mcost)
        scosts.append(scost)

    pylab.scatter(mcosts, scosts)
    pylab.savefig('pareto_mst_test.pdf', format='pdf')
    pylab.close()

if __name__ == '__main__':
    points = [(0, 0), (1, 1), (1, 1.1), (0, 0.1), (2, 2), (-1, -1), (-1, -1.1), (-1, 2), (-0.5, -0.5), (-0.5, 0.5), (0.5, 0.5), (1.1, 0.01)]
    root = (0, 0)
    pareto_plot(points, root)
