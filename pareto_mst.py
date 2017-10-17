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
from kruskal import *
from neuron_builder import build_neuron
from random_graphs import random_mst
from graph_utils import *

SKIP_TYPES = ['Unknown_neurotransmitter', 'Not_reported']

VIZ_TREES = False
MIN_NODES = 0
MAX_NODES = 3000
NO_CONTINUE = False

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
        pareto_tree = pareto_prim(H, alpha)
        viz_tree(pareto_tree, name + str(alpha), outdir=outdir)


    viz_tree(G, name + '_neural', outdir=outdir)

    sat_tree = satellite_tree(H)
    viz_tree(sat_tree, name + '_satellite', outdir=outdir)

    mst = min_spanning_tree(H)
    viz_tree(mst, name + '_mst', outdir=outdir)

def pareto_plot(G, name, cell_type, species, region, lab, outdir='figs',\
                output=True, viz_trees=VIZ_TREES, axon=False, compare=True):
    if not (nx.is_connected(G) and G.number_of_edges() == G.number_of_nodes() - 1):
        print "not a tree"
        return None

    assert G.number_of_nodes() > 0
   
    print "making graph"
    point_graph = None
    if NO_CONTINUE:
        point_graph = non_continue_subgraph(G)
    else:
        point_graph = complete_graph(G)
    print point_graph.number_of_nodes(), "points"
   
    sat_tree = satellite_tree(point_graph)
    span_tree = nx.minimum_spanning_tree(point_graph, weight='length')
    
    opt_scost = float(satellite_cost(sat_tree))
    normalize_scost = lambda cost : normalize_cost(cost, opt_scost)
    opt_mcost = float(mst_cost(span_tree))
    normalize_mcost = lambda cost : normalize_cost(cost, opt_mcost)
    
    mcosts1 = []
    scosts1 = []
    
    mcosts2 = []
    scosts2 = []

    mcosts3 = []
    scosts3 = []

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
            pareto_tree3 = sat_tree
        elif alpha == 1:
            pareto_tree1 = pareto_tree2 = span_tree
            pareto_tree3 = span_tree
        else:
            #print "Greedy"
            pareto_tree1 = pareto_prim(point_graph, alpha, axon=axon)
            #pareto_tree1 = pareto_prim_sandbox(point_graph, alpha, axon=axon)
            #print "khuller"
            if compare:
                beta = alpha_to_beta(alpha, opt_mcost, opt_scost)
                pareto_tree2 = khuller(point_graph, span_tree, sat_tree, beta)
            #print "genetic"
            #pareto_tree3 = pareto_genetic(point_graph, alpha)

        assert is_tree(pareto_tree1)
        if compare:
            assert is_tree(pareto_tree2)
        #assert is_tree(pareto_tree3)

        if (alpha != 0 and alpha != 1 and i % 5 == 0) and viz_trees:
            viz_tree(pareto_tree1, name + '-' + str(alpha), outdir=outdir)
        
        mcost1, scost1 = graph_costs(pareto_tree1)
        mcost1 = normalize_mcost(mcost1)
        scost1 = normalize_scost(scost1)
       
        if compare:
            mcost2, scost2 = graph_costs(pareto_tree2)
            mcost2 = normalize_mcost(mcost2)
            scost2 = normalize_scost(scost2)
        
        #mcost3, scost3 = graph_costs(pareto_tree3)
        #mcost3 = normalize_mcost(mcost3)
        #scost3 = normalize_scost(scost3)
        
        mcosts1.append(mcost1)
        scosts1.append(scost1)
       
        if compare:
            mcosts2.append(mcost2)
            scosts2.append(scost2)

        #mcosts3.append(mcost3)
        #scosts3.append(scost3)
        
        if alpha not in [0, 1] and compare:
            comparisons += 1
            if pareto_cost(mcost1, scost1, alpha) <= pareto_cost(mcost2, scost2, alpha):
                dominates += 1

    pylab.figure()
    pylab.plot(mcosts1, scosts1, c = 'b')
    pylab.scatter(mcosts1, scosts1, c='b', label='greedy mst')
    
    if compare:
        pylab.plot(mcosts2, scosts2, c = 'k')
        pylab.scatter(mcosts2, scosts2, c='k', label='khuller mst')
    
    #pylab.plot(mcosts3, scosts3, c = 'y')
    #pylab.scatter(mcosts3, scosts3, c='y', label='genetic')
    
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
        write_items = [name, cell_type, species, region, lab, point_graph.number_of_nodes()]
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

def pareto_plot_imaris(G, name, outdir='figs', viz_trees=VIZ_TREES, axon=False):
    if not is_tree(G):
        print "not a tree"
        return None

    assert G.number_of_nodes() > 0
   
    print "making graph"
    point_graph = None
    if NO_CONTINUE:
        point_graph = non_continue_subgraph(G)
    else:
        point_graph = complete_graph(G)
    print point_graph.number_of_nodes(), "points"
   
    sat_tree = satellite_tree(point_graph)
    span_tree = nx.minimum_spanning_tree(point_graph, weight='length')
    
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

    for i, alpha in enumerate(alphas):
        print "alpha", alpha
        pareto_tree1 = None
        pareto_tree2 = None
        if alpha == 0:
            pareto_tree1 = pareto_tree2 = sat_tree
            pareto_tree3 = sat_tree
        elif alpha == 1:
            pareto_tree1 = pareto_tree2 = span_tree
            pareto_tree3 = span_tree
        else:
            pareto_tree1 = pareto_prim(point_graph, alpha, axon=axon)
            beta = alpha_to_beta(alpha, opt_mcost, opt_scost)
            pareto_tree2 = khuller(point_graph, span_tree, sat_tree, beta)

        assert is_tree(pareto_tree1)
        assert is_tree(pareto_tree2)

        if (alpha != 0 and alpha != 1 and i % 5 == 0) and viz_trees:
            viz_tree(pareto_tree1, name + '-' + str(alpha), outdir=outdir)
        
        mcost1, scost1 = graph_costs(pareto_tree1)
        mcost1 = normalize_mcost(mcost1)
        scost1 = normalize_scost(scost1)
       
        mcost2, scost2 = graph_costs(pareto_tree2)
        mcost2 = normalize_mcost(mcost2)
        scost2 = normalize_scost(scost2)
         
        mcosts1.append(mcost1)
        scosts1.append(scost1)
       
        mcosts2.append(mcost2)
        scosts2.append(scost2) 

    pylab.figure()
    pylab.plot(mcosts1, scosts1, c = 'b')
    pylab.scatter(mcosts1, scosts1, c='b', label='prim')
    
    pylab.plot(mcosts2, scosts2, c = 'k')
    pylab.scatter(mcosts2, scosts2, c='k', label='khuller')
    
    pylab.plot(mcosts3, scosts3, c = 'r')
    pylab.scatter(mcosts3, scosts3, c='r', label='steiner')
    
    pylab.xlabel('spanning tree cost')
    pylab.ylabel('satellite cost')
    
    neural_mcost = normalize_mcost(mst_cost(G))
    neural_scost  = normalize_scost(satellite_cost(G))

    centroid_tree = centroid_mst(point_graph)

    centroid_mcost = normalize_mcost(mst_cost(centroid_tree))
    centroid_scost = normalize_scost(satellite_cost(centroid_tree))
    

    pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+', linewidths=15, label='centroid mst')

    #pylab.axhline(opt_scost)
    #pylab.axvline(opt_mcost)

    #pylab.legend()

    pylab.savefig('%s/pareto_front_%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

    viz_tree(G, name + str('_neural'), outdir=outdir) 
    viz_tree(sat_tree, name + str('_sat'), outdir=outdir)
    viz_tree(span_tree, name + str('_mst'), outdir=outdir)
    viz_tree(centroid_tree, name + '_centroid', outdir=outdir) 

def neuromorpho_plots(min_nodes=MIN_NODES, max_nodes=MAX_NODES):
    #directory = 'neuromorpho'
    directory = 'datasets'
    i = 0
    for cell_type in os.listdir(directory):
        if cell_type in SKIP_TYPES:
            continue
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
                            if not (min_nodes <= G.number_of_nodes() <= max_nodes):
                                continue
                            name = neuron[:-8] + str(i)
                            outdir = 'figs/%s/%s/%s/%s/%s' % (cell_type, species, region, lab, name)
                            outdir = outdir.replace(' ', '_')
                            os.system('mkdir -p %s' % outdir)
                            if len(os.listdir(outdir)) > 0:
                                continue
                            print species, lab, neuron
                            axon = i == 0
                            try:
                                pareto_plot(G, name, cell_type, species, region,\
                                            lab, outdir, axon=axon)
                            except RuntimeError as r:
                                print r
                                continue

def imaris_plots():
    for subdir in os.listdir('imaris'):
        if os.path.isdir('imaris/' + subdir):
            imfiles = []
            for fname in os.listdir('imaris/' + subdir):
                if 'Position' in fname:
                    imfiles.append('imaris/' + subdir + '/' + fname)
            G = read_imaris(imfiles, viz=False)
            outdir = 'imaris/' + subdir
            pareto_plot_imaris(G, subdir, outdir, viz_trees=True)


def neuron_builder_plots(rmin=0.5, rmax=1.5, rstep=0.01, num_iters=10):
    for i in xrange(num_iters):
        for radius in pylab.arange(rmax, rmin - rstep, -rstep):
            print "radius", radius
            G = build_neuron(radius)
            point_graph = complete_graph(G)
            print point_graph.number_of_nodes(), "points"

            print "sorting neighbors"
            sort_neighbors(point_graph)
           
            sat_tree = satellite_tree(point_graph)
            span_tree = nx.minimum_spanning_tree(point_graph, weight='length')
            
            opt_scost = float(satellite_cost(sat_tree))
            normalize_scost = lambda cost : normalize_cost(cost, opt_scost)
            opt_mcost = float(mst_cost(span_tree))
            normalize_mcost = lambda cost : normalize_cost(cost, opt_mcost)
            
            mcosts = []
            scosts = []
            
            delta = 0.01
            alphas = pylab.arange(0, 1 + delta, delta)
            for i, alpha in enumerate(alphas):
                print "alpha", alpha
                pareto_tree = None
                if alpha == 0:
                    pareto_tree = sat_tree
                elif alpha == 1:
                    pareto_tree = span_tree
                else:
                    pareto_tree = pareto_prim(point_graph, alpha)
    
                mcost, scost = graph_costs(pareto_tree)
                mcost = normalize_mcost(mcost)
                scost = normalize_scost(scost)
                mcosts.append(mcost)
                scosts.append(scost)

            neural_mcost = normalize_mcost(mst_cost(G))
            neural_scost  = normalize_scost(satellite_cost(G))
            
            pylab.figure()
            pylab.scatter(mcosts, scosts)
            pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x')
            pylab.plot(mcosts, scosts)
            fname = 'pareto_front%d' % len(os.listdir('neuron_builder'))
            pylab.savefig('%s/%s.pdf' % ('neuron_builder', fname), format='pdf')

            neural_dist, neural_index = pareto_dist(mcosts, scosts, neural_mcost, neural_scost)
            neural_closem = mcosts[neural_index] 
            neural_closes = scosts[neural_index]

            neural_alpha = alphas[neural_index]

            neural_cost = pareto_cost(neural_mcost, neural_scost, neural_alpha)
            neural_close_cost = pareto_cost(neural_closem, neural_closes, neural_alpha)
            if  neural_cost < neural_close_cost:
                neural_dist *= -1

            write_items = [radius, G.number_of_nodes(), neural_dist, neural_alpha]
            write_items = map(str, write_items)
            write_items = ', '.join(write_items)
            outfile = open('neuron_builder.csv', 'a')
            outfile.write('%s\n' % write_items)
            outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-min_nodes', type=int, default=MIN_NODES)
    parser.add_argument('-max_nodes', type=int, default=MAX_NODES)
    parser.add_argument('-a', '--algorithm', choices=['greedy', 'khuller'], default='greedy')
    parser.add_argument('-n', '--neuromorpho', action='store_true')
    parser.add_argument('-i', '--imaris', action='store_true')
    parser.add_argument('-b', '--neuron_builder', action='store_true')

    args = parser.parse_args()
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes
    algorithm = args.algorithm
    neuromorpho = args.neuromorpho
    imaris = args.imaris
    neuron_builder = args.neuron_builder

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
    if neuron_builder:
        neuron_builder_plots()

    #pareto_plot(thinstar_filename, 'thinstar', 'amacrine', 'rabbit', 'retina', 'miller', outdir='sandbox')

    #pareto_drawings(goldfish_filename, 'goldfish')
    #pareto_drawings(frog_filename, 'frog')
