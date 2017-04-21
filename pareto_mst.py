import numpy as np
import networkx as nx
from sys import argv
import matplotlib as mpl
mpl.use('agg')
import pylab

def get_neuron_points(filename, dim='3D'):
    #for arbor_type in ["2","3","4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
    H = nx.Graph()
    H_root = None
    HP2Coord = None
    for arbor_type in ["2", "3", "4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
        
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

        if G.number_of_nodes() > H.number_of_nodes():
            H = G
            HP2Coord = P2Coord
            H_root = root

                    #if UTIL.euclidean_dist(u,v,P2Coord) > 10:
                    #    print u,v
    return H, HP2Coord, H_root
    

def makes_cycle(u, v, node_to_forest):
    f1 = node_to_forest[u]
    f2 = node_to_forest[v]
    return f1 == f2
 
def combine_forests(u, v, node_to_forest, forest_to_nodes, forest_to_edges):
    f1 = node_to_forest[u]
    f2 = node_to_forest[v]
    if f1 == f2:
        return
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
    H = G.copy()
   
    pareto_mst.add_node(root)
    pareto_mst.node[root]['droot'] = 0
   
    graph_mcost = 0
    graph_scost = 0
 
    while pareto_mst.number_of_nodes() < H.number_of_nodes():
        best_edge = None
        best_cost = float("inf")

        #print "H edges", H.number_of_edges()

        candidate_edges = set()
        for u in pareto_mst.nodes():
            candidate_neighbors = H.neighbors(u)
            closest_neigbor = None
            closest_dist = float("inf")
            assert 'droot' in pareto_mst.node[u]
            for v in candidate_neighbors:
                if pareto_mst.has_node(v):
                    H.remove_edge(u, v)
                    #assert 'droot' in pareto_mst.node[v]
                elif G[u][v]['length'] < closest_dist:
                    closest_neighbor = v
                    closest_dist = H[u][v]['length']

            candidate_edges.add(tuple(sorted((u, closest_neighbor))))

        for i, (u, v) in enumerate(candidate_edges):
            #print "candidate edge", i, len(candidate_edges)

            scost = graph_scost           
            if pareto_mst.has_node(u):
                #assert 'droot' in pareto_mst.node[u]
                #assert not pareto_mst.has_node(v)
                scost += pareto_mst.node[u]['droot'] + G[u][v]['length']
            elif pareto_mst.has_node(v):
                #assert 'droot' in pareto_mst.node[v]
                #assert not pareto_mst.has_node(u)
                scost += pareto_mst.node[v]['droot'] + G[u][v]['length']
            else:
                raise ValueError('something is wrong')

            mcost = graph_mcost + G[u][v]['length']
 
            mcost *= alpha
            scost *= (1 - alpha)
            cost = mcost + scost
            
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
        if best_edge == None:
            break
        u, v = best_edge
        pareto_mst.add_edge(u, v)
        pareto_mst[u][v]['length'] = G[u][v]['length']
        
        #assert ('droot' in pareto_mst.node[u]) or ('droot' in pareto_mst.node[v]) 
        #assert ('droot' not in pareto_mst.node[u]) or ('droot' not in pareto_mst.node[v])

        if 'droot' in pareto_mst.node[u]:
            #assert 'droot' not in pareto_mst.node[v]
            pareto_mst.node[v]['droot'] = pareto_mst.node[u]['droot'] + pareto_mst[u][v]['length']
        elif 'droot' in pareto_mst.node[v]:
            #assert 'droot' not in pareto_mst.node[u]
            pareto_mst.node[u]['droot'] = pareto_mst.node[v]['droot'] + pareto_mst[u][v]['length']
        else:
            raise ValueError('something is really wrong')


        H.remove_edge(u, v)
    
    return pareto_mst 

'''
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
        print "number of candidate edge", len(candidate_edges)
        best_edge = None
        best_cost = float('inf')
        for i, (u, v) in enumerate(candidate_edges):
            print "candidate edge", i
            u, v = sorted((u, v))
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
'''

def point_dist(p1, p2):
    assert len(p1) == len(p2)
    sq_dist = 0
    for i in xrange(len(p1)):
        x1, x2 = p1[i], p2[i]
        sq_dist += (x1 - x2) ** 2
    return sq_dist ** 0.5

def pareto_mst(points, root, alpha):
    G = nx.Graph()
    print "making graph"
    for i in xrange(len(points)):
        for j in xrange(len(points)):
            if i == j:
                continue
            p1, p2 = points[i], points[j]
            G.add_edge(p1, p2)
            length = point_dist(p1, p2)
            G[p1][p2]['length'] = length
    print "computing pareto mst"
    mst = pareto_kruskal(G, root, alpha)
    return mst

def centroid_mst(points, root):
    centroid = np.zeros(len(root))
    for point in points:
        assert len(point) == len(root)
        if point != root:
            centroid += point
    centroid /= len(points) - 1

    root_cdist = point_dist(root, centroid)

    mcost = root_cdist
    scost = 0
    
    for point in points:
        if point != root:
            cdist = point_dist(point, centroid)
            scost += cdist
            mcost += cdist + root_cdist

    return mcost, scost



def pareto_plot(filename, name):
    G, P2Coord, root_node = get_neuron_points(filename)
    points = P2Coord.values()
    root = P2Coord[root_node]
    
    mcosts = []
    scosts = []
    delta = 0.01
    for alpha in np.arange(0, 1 + delta, delta):
        print "alpha", alpha
        mst = pareto_mst(points, root, alpha)
        mcost = mst_cost(mst)
        scost = satellite_cost(mst, root)
        
        mcosts.append(mcost)
        scosts.append(scost)

    pylab.plot(mcosts, scosts, c = 'b')
    pylab.scatter(mcosts, scosts, c='b')
    pylab.xlabel('spanning tree cost')
    pylab.ylabel('satellite cost')
    neural_mst = nx.Graph()
    for u, v in G.edges_iter():
        p1, p2 = P2Coord[u], P2Coord[v]
        neural_mst.add_edge(p1, p2)
        neural_mst[p1][p2]['length'] = point_dist(p1, p2)

    neural_length = mst_cost(neural_mst)
    neural_droot  = satellite_cost(neural_mst, root)
    
    pylab.scatter([neural_length], [neural_droot], c='r', marker='x', linewidths=15)

    centroid_droot, centroid_length = centroid_mst(points, root)
    pylab.scatter([centroid_length], [centroid_droot], c='g', marker='+', linewidths=15)

    xmin = ymin = 0
    xmax = max(mcosts + [neural_length])
    ymax = max(scosts + [neural_droot])


    curr_ax = pylab.gca()
    #curr_ax.set_xlim([xmin,xmax + 100])
    #curr_ax.set_ylim([ymin,ymax + 100])

    pylab.savefig('pareto_mst_%s.pdf' % name, format='pdf')
    pylab.close()

if __name__ == '__main__':
    points = [(0, 0), (1, 1), (1, 1.1), (0, 0.1), (2, 2), (-1, -1), (-1, -1.1), (-1, 2), (-0.5, -0.5), (-0.5, 0.5), (0.5, 0.5), (1.1, 0.01)]
    root = (0, 0)
    #pareto_plot(points, root, 'test')

    #filename = argv[1]
    frog_filename = 'neuromorpho/frog/birinyi/GEN1.CNG.swc'
    goldfish_filename= 'neuromorpho/goldfish/stevens/G4-19g-1.CNG.swc'
    pig_filename = 'neuromorpho/pig/johnson/Pig288-DG-In-Neuron1.CNG.swc'
    agouti_filename = 'neuromorpho/agouti/guerra da rocha/cco6lam06cel05pa.CNG.swc'
    celegans_filename = 'neuromorpho/celegans/openworm/SABVR.CNG.swc' 
    
    #points = P2Coord.values()

    #print P2Coord.keys()

    #root_point = P2Coord[root]

    pareto_plot(goldfish_filename, 'goldfish')
    pareto_plot(pig_filename, 'pig')
    pareto_plot(agouti_filename, 'agouti')
    pareto_plot(celegans_filename, 'celegans')
    pareto_plot(frog_filename, 'frog')
    
    '''
    pareto_plot(G, points, root_point, P2Coord, 'pig')
    pareto_plot(G, points, root_point, P2Coord, 'pig')
    pareto_plot(G, points, root_point, P2Coord, 'pig')
    pareto_plot(G, points, root_point, P2Coord, 'pig')
    pareto_plot(G, points, root_point, P2Coord, 'pig')
    '''
