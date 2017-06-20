import numpy as np
import networkx as nx
from sys import argv
import matplotlib as mpl
mpl.use('agg')
import pylab
import os
from random import shuffle
from itertools import combinations
import argparse
from khuller import khuller
from cost_functions import *
from neuron_utils import *
from read_imaris import *
from pareto_functions import *

VIZ_TREES = False
MIN_NODES = 0
MAX_NODES = 3000
NO_CONTINUE = False

def get_neuron_points(filename, dim='3D'):
    #for arbor_type in ["2","3","4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
    graphs = []
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
                    coord = None
                    if dim == "3D":
                        coord = (float(cols[2]),float(cols[3]),float(cols[4]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        coord = (float(cols[2]),float(cols[3]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]))
                    else:
                        assert False
                    G.node[root]['coord'] = coord
                    G.graph['root'] = root


                else:
                    u,v = int(cols[6]),int(cols[0])
                    assert u in G #and u in P2Coord
                    assert not G.has_edge(u,v)
            
                    coord = None
                    if dim == "3D":
                        coord = (float(cols[2]),float(cols[3]),float(cols[4]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        coord = (float(cols[2]),float(cols[3]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]))
                    else:
                        assert False
                    G.add_edge(u, v)
                    
                    G.node[v]['coord'] = coord
                    
                    G[u][v]['length'] = point_dist(G.node[u]['coord'], G.node[v]['coord'])

        assert 'root' in G.graph
        label_points(G)
        graphs.append(G)

    return graphs

def initialize_lengths(G):
    for u, v in G.edges_iter():
        p1, p2 = G.node[u]['coord'], G.node[v]['coord']
        G[u][v]['length'] = point_dist(p1, p2)

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
    mst = G.copy()
    mst.remove_edges_from(G.edges())
    for u, v in mst_edges:
        mst.add_edge(u, v)
        mst[u][v]['length'] = G[u][v]['length']
    return mst

def centroid_mst(G):
    cent_mst = G.copy()
    cent_mst.remove_edges_from(G.edges())
    
    centroidp = centroid(G)
    cent_mst.add_node('centroid')
    cent_mst.node['centroid']['label'] = 'centroid'
    cent_mst.node['centroid']['coord'] = centroidp
    for u in G.nodes_iter():
        cent_mst.add_edge(u, 'centroid')
        cent_mst[u]['centroid']['length'] = point_dist(cent_mst.node[u]['coord'], centroidp)
    return cent_mst

def centroid_mst_costs(G):
    root = G.graph['root']
    centroidp = centroid(G)
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
    for p1, p2 in combinations(points, 2):
        point_graph.add_edge(p1, p2)
        length = point_dist(p1, p2)
        point_graph[p1][p2]['length'] = length
    return point_graph

def complete_graph(G):
    H = G.copy()
    H.remove_edges_from(H.edges())
    for u, v in combinations(H.nodes(), 2):
        H.add_edge(u, v)
        H[u][v]['length'] = point_dist(H.node[u]['coord'], H.node[v]['coord'])
    return H

def pareto_dist(pareto_mcosts, pareto_scosts, mcost, scost):
    best_dist = float("inf")
    best_index = None

    assert len(pareto_mcosts) == len(pareto_scosts)

    for i in xrange(len(pareto_mcosts)):
        pareto_mcost = pareto_mcosts[i]
        pareto_scost = pareto_scosts[i]

        dist = point_dist((pareto_mcost, pareto_scost), (mcost, scost))
        if dist < best_dist:
            best_dist = dist
            best_index = i

    return best_dist, best_index

def satellite_tree(G):
    root = G.graph['root']
    
    satellite = G.copy()
    satellite.remove_edges_from(G.edges())

    for u in satellite.nodes_iter():
        if u != root:
            satellite.add_edge(u, root)
            p1, p2 = satellite.node[u]['coord'], satellite.node[root]['coord']
            satellite[u][root]['length'] = point_dist(p1, p2)
    
    return satellite

def best_satellite_cost(G):
    best_cost = 0
    root = G.graph['root']
    root_coord = G.node[root]['coord']
    for u in G.nodes_iter():
        best_cost += point_dist(root_coord, G.node[u]['coord'])
    return best_cost

def min_spanning_tree(G):
    sorted_edges = sorted(G.edges(), key= lambda x : G[x[0]][x[1]]['length'])
    mst_edges = kruskal(G.nodes(), sorted_edges)
    mst = G.copy()
    mst.remove_edges_from(mst.edges())
    for u, v in mst_edges:
        mst.add_edge(u, v)
        mst[u][v]['length'] = G[u][v]['length']
    return mst

def best_mst_cost(G):
    sorted_edges = sorted(G.edges(), key= lambda x : G[x[0]][x[1]]['length'])
    mst_edges = kruskal(G.nodes(), sorted_edges)
    best_cost = 0
    for u, v in mst_edges:
        best_cost += G[u][v]['length']
    return best_cost 

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

def pareto_drawings(filename, name, outdir='drawings'): 
    G = get_neuron_points(filename)
    H = non_continue_subgraph(G)
    if H.number_of_nodes() <= 0 or H.number_of_nodes() > 1000:
        return None
 
    outdir = '%s/%s' % (outdir, name)
    os.system('mkdir -p %s' % outdir)

    delta = 0.01
    alphas = np.arange(0, 1 + delta, delta)
    for alpha in alphas:
        print "alpha", alpha
        pareto_tree = pareto_kruskal(H, alpha)
        viz_tree(pareto_tree, name + str(alpha), outdir=outdir)


    viz_tree(G, name + '_neural', outdir=outdir)

    sat_tree = satellite_tree(H)
    viz_tree(sat_tree, name + '_satellite', outdir=outdir)

    mst = min_spanning_tree(H)
    viz_tree(mst, name + '_mst', outdir=outdir)

def normalize_cost(cost, opt_cost):
    return 1 - (opt_cost / cost)

def sort_neighbors(G):
    for u in G.nodes_iter():
        G.node[u]['close_neighbors'] = sorted(G.neighbors(u), key = lambda v : G[u][v]['length'])

def pareto_plot(G, name, cell_type, species, region, lab, outdir='figs',\
                min_nodes=MIN_NODES, max_nodes=MAX_NODES, output=True, viz_trees=VIZ_TREES):
    if not (nx.is_connected(G) and G.number_of_edges() == G.number_of_nodes() - 1):
        print "not a tree"
        return None

    assert G.number_of_nodes() > 0
    print G.number_of_nodes(), "nodes"
    if G.number_of_nodes() <= min_nodes or G.number_of_nodes() > max_nodes:
        return None
   
    print "making graph"
    point_graph = None
    if NO_CONTINUE:
        point_graph = non_continue_subgraph(G)
    else:
        point_graph = complete_graph(G)
    print point_graph.number_of_nodes(), "points"
   
    sat_tree = satellite_tree(point_graph)
    span_tree = min_spanning_tree(point_graph)
    
    opt_scost = float(satellite_cost(sat_tree))
    normalize_scost = lambda cost : normalize_cost(cost, opt_scost)
    opt_mcost = float(mst_cost(span_tree))
    normalize_mcost = lambda cost : normalize_cost(cost, opt_mcost)
    
    mcosts1 = []
    scosts1 = []
    mcosts2 = []
    scosts2 = []

    delta = 0.01
    alphas = np.arange(0, 1 + delta, delta)
    print "sorting neighbors"
    sort_neighbors(point_graph)

    comparisons = 0
    dominates = 0
    for i, alpha in enumerate(alphas):
        print "alpha", alpha
        pareto_tree1 = None
        pareto_tree2 = None
        if alpha == 0:
            pareto_tree1 = pareto_tree2 = sat_tree
        elif alpha == 1:
            pareto_tree1 = pareto_tree2 = span_tree
        else:
            pareto_tree1 = pareto_kruskal(point_graph, alpha)
            pareto_tree2 = khuller(point_graph, span_tree, sat_tree, 1.0 / (1 - alpha))

        assert is_tree(pareto_tree1)
        assert is_tree(pareto_tree2)

        if (alpha != 0 and alpha != 1 and i % 5 == 0) and viz_trees:
            viz_tree(pareto_tree1, name + '-' + str(alpha), outdir=outdir)
        
        mcost1 = normalize_mcost(mst_cost(pareto_tree1))
        scost1 = normalize_scost(satellite_cost(pareto_tree1))
        
        mcost2 = normalize_mcost(mst_cost(pareto_tree2))
        scost2 = normalize_scost(satellite_cost(pareto_tree2))
        
        mcosts1.append(mcost1)
        scosts1.append(scost1)
        
        mcosts2.append(mcost2)
        scosts2.append(scost2)

        if alpha not in [0, 1]:
            comparisons += 1
            if pareto_cost(mcost1, scost1, alpha) <= pareto_cost(mcost2, scost2, alpha):
                dominates += 1

    pylab.figure()
    pylab.plot(mcosts1, scosts1, c = 'b')
    pylab.scatter(mcosts1, scosts1, c='b', label='greedy mst')
    
    pylab.plot(mcosts2, scosts2, c = 'k')
    pylab.scatter(mcosts2, scosts2, c='k', label='khuller mst')
    
    pylab.xlabel('spanning tree cost')
    pylab.ylabel('satellite cost')
    
    neural_mcost = normalize_mcost(mst_cost(G))
    neural_scost  = normalize_scost(satellite_cost(G))

    neural_dist, neural_index = pareto_dist(mcosts1, scosts1, neural_mcost, neural_scost)
    neural_closem = mcosts1[neural_index] 
    neural_closes = scosts1[neural_index]
    neural_alpha = alphas[neural_index]

    pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x', linewidths=15, label='neural mst')
    #pylab.plot([neural_mcost, neural_closem], [neural_scost, neural_closes], c='r', linestyle='--')

    centroid_tree = centroid_mst(point_graph)

    centroid_mcost = normalize_mcost(mst_cost(centroid_tree))
    centroid_scost = normalize_scost(satellite_cost(centroid_tree))
    
    centroid_dist, centroid_index = pareto_dist(mcosts1, scosts1, centroid_mcost, centroid_scost)
    centroid_closem = mcosts1[centroid_index]
    centroid_closes = scosts1[centroid_index]
    centroid_alpha = alphas[centroid_index]

    pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+', linewidths=15, label='centroid mst')
    #pylab.plot([centroid_mcost, centroid_closem], [centroid_scost, centroid_closes], c='g', linestyle='--')

    #pylab.axhline(opt_scost)
    #pylab.axvline(opt_mcost)

    #pylab.legend()

    ntrials = 20
    successes = 0
    total_rand_dist = 0.0
    rand_mcosts = []
    rand_scosts = []
    for i in xrange(ntrials):
        rand_mst = random_mst(point_graph)
        
        rand_mcost = normalize_mcost(mst_cost(rand_mst))
        rand_scost = normalize_scost(satellite_cost(rand_mst))
        rand_mcosts.append(rand_mcost)
        rand_scosts.append(rand_scost)
        
        rand_dist, rand_index = pareto_dist(mcosts1, scosts1, rand_mcost, rand_scost)
        total_rand_dist += rand_dist
        rand_closem = mcosts1[rand_index]
        rand_closes = scosts1[rand_index]
        rand_alpha = alphas[rand_index]
        if rand_dist < neural_dist:
            successes += 1

    pylab.scatter(rand_mcosts, rand_scosts, c='m', marker='o', label='random mst')

    pylab.savefig('%s/pareto_front_%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

    viz_tree(G, name + str('_neural'), outdir=outdir) 
    viz_tree(sat_tree, name + str('_sat'), outdir=outdir)
    viz_tree(span_tree, name + str('_mst'), outdir=outdir)
    viz_tree(centroid_tree, name + '_centroid', outdir=outdir) 

    mean_rand_dist = total_rand_dist / ntrials
   
    f = open('pareto_mst.csv', 'a')
    if output:
        write_items = [name, cell_type, species, region, lab]
        write_items = map(str, write_items)
        write_items.append(neural_alpha)
        write_items += [neural_dist, centroid_dist, mean_rand_dist]
        write_items += [ntrials, successes]
        write_items += [comparisons, dominates]
        write_items = map(str, write_items)
        write_items = ', '.join(write_items)
        f.write('%s\n' % write_items)

    #pylab.scatter([rand_mcost], [rand_scost], c='m', marker='*', linewidths=15)
    #pylab.plot([rand_mcost, rand_closem], [rand_scost, rand_closes], c='m', linestyle='-')

def neuromorpho_plots(min_nodes=MIN_NODES, max_nodes=MAX_NODES):
    #directory = 'neuromorpho'
    directory = 'datasets'
    i = 0
    for cell_type in os.listdir(directory):
        for species in os.listdir(directory + '/' + cell_type):
            for region in os.listdir(directory + '/' + cell_type + '/' + species):
                for lab in os.listdir(directory + "/" + cell_type + '/' + species+ '/' + region):
                    for neuron in os.listdir(directory + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = directory + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron
                        
                        if neuron[-8:] != ".CNG.swc": 
                            continue
                        
                        try:
                            graphs = get_neuron_points(filename)
                        except AssertionError:
                            continue

                        for i, G in enumerate(graphs):
                            if G.number_of_nodes() <= 1:
                                continue
                            name = neuron[:-8] + str(i)
                            outdir = 'figs/%s/%s/%s/%s/%s' % (cell_type, species, region, lab, name)
                            outdir = outdir.replace(' ', '_')
                            os.system('mkdir -p %s' % outdir)
                            if len(os.listdir(outdir)) > 0:
                                continue
                            print species, lab, neuron
                            pareto_plot(G, name, cell_type, species, region,\
                                        lab, outdir, min_nodes, max_nodes)

def imaris_plots():
    for subdir in os.listdir('imaris'):
        if os.path.isdir('imaris/' + subdir):
            imfiles = []
            for fname in os.listdir('imaris/' + subdir):
                if 'Position' in fname:
                    imfiles.append('imaris/' + subdir + '/' + fname)
            G = read_imaris(imfiles, viz=False)
            outdir = 'imaris/' + subdir
            pareto_plot(G, subdir, None, None, None, None, outdir, output=False, viz_trees=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-min_nodes', type=int, default=MIN_NODES)
    parser.add_argument('-max_nodes', type=int, default=MAX_NODES)
    parser.add_argument('-a', '--algorithm', choices=['greedy', 'khuller'], default='greedy')
    parser.add_argument('-n', '--neuromorpho', action='store_true', default=False)
    parser.add_argument('-i', '--imaris', action='store_true', default=False)

    args = parser.parse_args()
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes
    algorithm = args.algorithm
    neuromorpho = args.neuromorpho
    imaris = args.imaris

    points = [(0, 0), (1, 1), (1, 1.1), (0, 0.1), (2, 2), (-1, -1), (-1, -1.1), (-1, 2), (-0.5, -0.5), (-0.5, 0.5), (0.5, 0.5), (1.1, 0.01)]
    root = (0, 0)
    #pareto_plot(points, root, 'test')

    #filename = argv[1]
    frog_filename = 'neuromorpho/frog/birinyi/GEN1.CNG.swc'
    goldfish_filename= 'neuromorpho/goldfish/stevens/G4-19g-1.CNG.swc'
    pig_filename = 'neuromorpho/pig/johnson/Pig288-DG-In-Neuron1.CNG.swc'
    agouti_filename = 'neuromorpho/agouti/guerra da rocha/cco6lam06cel05pa.CNG.swc'
    celegans_filename = 'neuromorpho/celegans/openworm/SABVR.CNG.swc' 
    thinstar_filename = 'datasets/amacrine/rabbit/retina/miller/THINSTAR.CNG.swc'
    
    #points = P2Coord.values()

    #print P2Coord.keys()

    #root_point = P2Coord[root]

    #pareto_plot(goldfish_filename, 'goldfish', 'cell_type1', 'species1',\
    #                               'region1', 'lab1', 'figs/',\
    #                               0, 1000, 'khuller')
    
    #pareto_plot(pig_filename, 'pig')
    #pareto_plot(agouti_filename, 'agouti')
    #pareto_plot(celegans_filename, 'celegans')
    #pareto_plot(frog_filename, 'frog')
    
    if imaris:
        imaris_plots()
    if neuromorpho:
        neuromorpho_plots(min_nodes, max_nodes)

    #pareto_plot(thinstar_filename, 'thinstar', 'amacrine', 'rabbit', 'retina', 'miller', outdir='sandbox')

    #pareto_drawings(goldfish_filename, 'goldfish')
    #pareto_drawings(frog_filename, 'frog')
