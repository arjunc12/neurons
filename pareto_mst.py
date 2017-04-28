import numpy as np
import networkx as nx
from sys import argv
import matplotlib as mpl
mpl.use('agg')
import pylab
import os
from random import shuffle

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

def random_mst(G):
    rand_edges = G.edges()
    shuffle(rand_edges)
    mst_edges = kruskal(G.nodes(), rand_edges)
    mst = nx.Graph()
    for u, v in mst_edges:
        mst.add_edge(u, v)
        mst[u][v]['length'] = point_dist(u, v)
    return mst

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

    print "sorting neighbors"
    closest_neighbors = {}
    for u in G.nodes():
        closest_neighbors[u] = sorted(G.neighbors(u), key=lambda v : G[u][v]['length'])

    while pareto_mst.number_of_nodes() < H.number_of_nodes():
        best_edge = None
        best_cost = float("inf")

        candidate_edges = set()
        for u in pareto_mst.nodes():
            assert 'droot' in pareto_mst.node[u]
             
            invalid_neighbors = []
            closest_neighbor = None
            for i in xrange(len(closest_neighbors[u])):
                v = closest_neighbors[u][i]
                if pareto_mst.has_node(v):
                    invalid_neighbors.append(v)
                else:
                    closest_neighbor = v
                    break

            for n in invalid_neighbors:
                closest_neighbors[u].remove(n)

            if closest_neighbor != None:
                candidate_edges.add(tuple(sorted((u, closest_neighbor))))

        for i, (u, v) in enumerate(candidate_edges):

            scost = graph_scost           
            if pareto_mst.has_node(u):
                scost += pareto_mst.node[u]['droot'] + G[u][v]['length']
            elif pareto_mst.has_node(v):
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
        

        if 'droot' in pareto_mst.node[u]:
            pareto_mst.node[v]['droot'] = pareto_mst.node[u]['droot'] + pareto_mst[u][v]['length']
        elif 'droot' in pareto_mst.node[v]:
            pareto_mst.node[u]['droot'] = pareto_mst.node[v]['droot'] + pareto_mst[u][v]['length']
        else:
            raise ValueError('something is really wrong')


        H.remove_edge(u, v)
        closest_neighbors[u].remove(v)
        closest_neighbors[v].append(u)
    
    return pareto_mst 

def point_dist(p1, p2):
    assert len(p1) == len(p2)
    sq_dist = 0
    for i in xrange(len(p1)):
        x1, x2 = p1[i], p2[i]
        sq_dist += (x1 - x2) ** 2
    return sq_dist ** 0.5

def centroid(points, root=None):
    centroid = np.zeros(len(root))
    for point in points:
        assert len(point) == len(root)
        if point != root:
            centroid += point
    centroid /= len(points) - 1
    return centroid

def centroid_mst(points, root):
    G = nx.Graph()
    centroidp = centroid(points, root)
    for point in points:
        G.add_edge(point, centroid)
        G[point][centroid]['length'] = point_dist(point, centroidp)
    return G

def centroid_mst_costs(points, root):
    centroidp = centroid(points, root)
    root_cdist = point_dist(root, centroidp)

    mcost = root_cdist
    scost = 0
    
    for point in points:
        if point != root:
            cdist = point_dist(point, centroidp)
            mcost += cdist
            scost += cdist + root_cdist

    return mcost, scost

def points_to_graph(points):
    point_graph = nx.Graph()
    for i in xrange(len(points)):
        for j in xrange(len(points)):
            if i == j:
                continue
            p1, p2 = points[i], points[j]
            point_graph.add_edge(p1, p2)
            length = point_dist(p1, p2)
            point_graph[p1][p2]['length'] = length
    return point_graph
   
def pareto_dist(pareto_mcosts, pareto_scosts, mcost, scost):
    best_pareto_mcost = None
    best_pareto_scost = None
    best_dist = float("inf")

    assert len(pareto_mcosts) == len(pareto_scosts)

    for i in xrange(len(pareto_mcosts)):
        pareto_mcost = pareto_mcosts[i]
        pareto_scost = pareto_scosts[i]

        dist = point_dist((pareto_mcost, pareto_scost), (mcost, scost))
        if dist < best_dist:
            best_dist = dist
            best_pareto_mcost = pareto_mcost
            best_pareto_scost = pareto_scost

    return best_pareto_mcost, best_pareto_scost, best_dist

def best_satellite_cost(root, points):
    best_cost = 0
    for point in points:
        best_cost += point_dist(root, point)
    return best_cost

def best_mst_cost(G):
    sorted_edges = sorted(G.edges(), key= lambda x : G[x[0]][x[1]]['length'])
    mst_edges = kruskal(G.nodes(), sorted_edges)
    best_cost = 0
    for u, v in mst_edges:
        best_cost += G[u][v]['length']
    return best_cost

def pareto_plot(filename, name, outdir=None):
    G, P2Coord, root_node = get_neuron_points(filename)
    points = P2Coord.values()
    root = P2Coord[root_node]
    print G.number_of_nodes(), "nodes"
    if G.number_of_nodes() <= 1000 or G.number_of_nodes() > 2000:
        return None
   
    print "making graph"
    point_graph = points_to_graph(points)
    
    mcosts = []
    scosts = []
    delta = 0.01
    for alpha in np.arange(0, 1 + delta, delta):
        print "alpha", alpha
        mst = pareto_kruskal(point_graph, root, alpha)
        mcost = mst_cost(mst)
        scost = satellite_cost(mst, root)
        
        mcosts.append(mcost)
        scosts.append(scost)

    pylab.plot(mcosts, scosts, c = 'b')
    pylab.scatter(mcosts, scosts, c='b', label='pareto mst')
    pylab.xlabel('spanning tree cost')
    pylab.ylabel('satellite cost')
    neural_mst = nx.Graph()
    for u, v in G.edges_iter():
        p1, p2 = P2Coord[u], P2Coord[v]
        neural_mst.add_edge(p1, p2)
        neural_mst[p1][p2]['length'] = point_dist(p1, p2)

    neural_mcost = mst_cost(neural_mst)
    neural_scost  = satellite_cost(neural_mst, root)

    neural_closem, neural_closes, neural_dist = pareto_dist(mcosts, scosts, neural_mcost, neural_scost)

    pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x', linewidths=15, label='neural mst')
    pylab.plot([neural_mcost, neural_closem], [neural_scost, neural_closes], c='r', linestyle='--')

    centroid_mcost, centroid_scost = centroid_mst_costs(points, root)
    centroid_closem, centroid_closes, centroid_dist = pareto_dist(mcosts, scosts, centroid_mcost, centroid_scost)

    pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+', linewidths=15, label='centroid mst')
    pylab.plot([centroid_mcost, centroid_closem], [centroid_scost, centroid_closes], c='g', linestyle='--')

    opt_scost = best_satellite_cost(root, points)
    opt_mcost = best_mst_cost(point_graph)

    #pylab.axhline(opt_scost)
    #pylab.axvline(opt_mcost)

    if outdir == None:
        outdir = 'figs'

    pylab.savefig('%s/pareto_mst_%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

    f = open('pareto_mst.csv', 'a')
    for i in xrange(100):
        rand_mst = random_mst(point_graph)
        rand_mcost = mst_cost(rand_mst)
        rand_scost = satellite_cost(rand_mst, root)
        rand_closem, rand_closes, rand_dist = pareto_dist(mcosts, scosts, rand_mcost, rand_scost)
        f.write('%s, %d, %d\n' % (name, neural_dist, rand_dist))



    #pylab.scatter([rand_mcost], [rand_scost], c='m', marker='*', linewidths=15)
    #pylab.plot([rand_mcost, rand_closem], [rand_scost, rand_closes], c='m', linestyle='-')


def neuromorpho_plots(plot_species=None):
    directory = 'neuromorpho'
    i = 0
    for species in os.listdir(directory):
        if plot_species != None and species not in plot_species:
            continue
        for lab in os.listdir(directory + "/" + species):
            for neuron in os.listdir(directory + "/" + species + "/" + lab):
                filename = directory + "/" + species + "/" + lab + "/" + neuron
                if neuron[-4:] != ".swc": 
                    continue
                name = neuron[:-4]
                outdir = 'figs/%s/%s' % (species, lab)
                outdir = outdir.replace(' ', '_')
                os.system('mkdir -p %s' % outdir)

                print species, lab, neuron

                pareto_plot(filename, neuron, outdir)


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

    #pareto_plot(goldfish_filename, 'goldfish')
    #pareto_plot(pig_filename, 'pig')
    #pareto_plot(agouti_filename, 'agouti')
    #pareto_plot(celegans_filename, 'celegans')
    #pareto_plot(frog_filename, 'frog')
    neuromorpho_plots()
