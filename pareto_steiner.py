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

SKIP_TYPES = ['Unknown_neurotransmitter', 'Not_reported']

VIZ_TREES = False
MIN_NODES = 0
MAX_NODES = 3000
NO_CONTINUE = False

def pareto_plot_neuromorpho(G, name, cell_type, species, region, lab,\
                            outdir='steiner_figs', output=True,\
                            viz_trees=VIZ_TREES, axon=False):
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
    opt_scost = satellite_cost(sat_tree)
    span_tree = nx.minimum_spanning_tree(point_graph, weight='length')

    '''
    opt_scost = float(satellite_cost(sat_tree))
    normalize_scost = lambda cost : normalize_cost(cost, opt_scost)
    opt_mcost = float(mst_cost(span_tree))
    normalize_mcost = lambda cost : normalize_cost(cost, opt_mcost)
    '''
    
    mcosts = []
    scosts = []
    costs = []

    delta = 0.01
    alphas = np.arange(delta, 1, delta)
    print "sorting neighbors"
    sort_neighbors(point_graph)

    comparisons = 0
    dominates = 0

    os.system("mkdir -p %s" % outdir)

    for i, alpha in enumerate(alphas):
        print "alpha", alpha
        pareto_tree = pareto_steiner(point_graph, alpha, axon=axon)
            
        assert is_tree(pareto_tree)
        mcost, scost = graph_costs(pareto_tree)
        cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
        if scost < opt_scost:
            print scost, opt_scost, opt_scost - scost
            print "dist error", dist_error(pareto_tree)
            print "scost error", scost_error(pareto_tree)
        assert scost >= opt_scost
        #check_dists(pareto_tree)

        if (i % 5 == 0) and viz_trees:
            viz_tree(pareto_tree1, name + '-' + str(alpha), outdir=outdir)
        
        mcosts.append(mcost)
        scosts.append(scost)
        costs.append(cost)

    min_mcost = min(mcosts)
    max_mcost = max(mcosts)
    min_scost = opt_scost
    max_scost = max(scosts)

    normalize_scost = make_normalize_function(min_scost)
    normalize_mcost = make_normalize_function(min_mcost)

    #normalized_mcosts = normalize_mcost(pylab.array(mcosts))
    normalized_mcosts = map(normalize_mcost, mcosts)
    #normalized_scosts = normalize_scost(pylab.array(scosts))
    normalized_scosts = map(normalize_scost, scosts)

    pylab.figure()
    pylab.plot(mcosts, scosts, c = 'b')
    pylab.scatter(mcosts, scosts, c='b', label='greedy steiner')
     
    pylab.xlabel('steiner tree cost')
    pylab.ylabel('satellite cost')
    
    #neural_mcost = normalize_mcost(mst_cost(G))
    #neural_scost  = normalize_scost(satellite_cost(G))
    neural_mcost, neural_scost = graph_costs(G)

    neural_dist, neural_index = pareto_dist(mcosts, scosts, neural_mcost,\
                                            neural_scost)
    neural_closem = mcosts[neural_index] 
    neural_closes = scosts[neural_index]
    neural_alpha = alphas[neural_index]

    neural_dist2, neural_index2 = pareto_dist(normalized_mcosts,\
                                              normalized_scosts, neural_mcost,\
                                              neural_scost)
    neural_alpha2 = alphas[neural_index]

    assert neural_alpha == neural_alpha2

    pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x', linewidths=15, label='neural')

    centroid_tree = centroid_mst(point_graph)

    #centroid_mcost = normalize_mcost(mst_cost(centroid_tree))
    #centroid_scost = normalize_scost(satellite_cost(centroid_tree))
    centroid_mcost, centroid_scost = graph_costs(centroid_tree)
    
    centroid_dist, centroid_index = pareto_dist(mcosts, scosts, centroid_mcost, centroid_scost)
    centroid_closem = mcosts[centroid_index]
    centroid_closes = scosts[centroid_index]
    centroid_alpha = alphas[centroid_index]

    pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+', linewidths=15, label='centroid mst')
    #pylab.plot([centroid_mcost, centroid_closem], [centroid_scost, centroid_closes], c='g', linestyle='--')

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
        
        rand_dist, rand_index = pareto_dist(mcosts, scosts, rand_mcost, rand_scost)
        total_rand_dist += rand_dist
        rand_closem = mcosts[rand_index]
        rand_closes = scosts[rand_index]
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
   
    f = open('pareto_steiner.csv', 'a')
    if output:
        write_items = [name, cell_type, species, region, lab, point_graph.number_of_nodes()]
        write_items.append(neural_alpha)
        write_items += [neural_dist, centroid_dist, mean_rand_dist]
        write_items += [ntrials, successes]
        write_items = map(str, write_items)
        write_items = ', '.join(write_items)
        f.write('%s\n' % write_items)

    f.close()

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
    
    mcosts1 = []
    scosts1 = []
    
    mcosts2 = []
    scosts2 = []

    delta = 0.01
    alphas = np.arange(delta, 1, delta)
    print "sorting neighbors"
    sort_neighbors(point_graph)

    os.system("mkdir -p %s" % outdir)

    sat_tree = satellite_tree(point_graph)
    opt_scost = float(satellite_cost(sat_tree))
    span_tree = nx.minimum_spanning_tree(point_graph, weight='length')

    for i, alpha in enumerate(alphas):
        print "alpha", alpha
        pareto_tree1 = pareto_steiner(point_graph, alpha, axon=axon)
        assert is_tree(pareto_tree1)
        mcost1, scost1 = graph_costs(pareto_tree1)
        assert scost1 >= opt_scost


        if (i % 5 == 0) and viz_trees:
            viz_tree(pareto_tree1, name + '-' + str(alpha), outdir=outdir)
        
        pareto_tree2 = pareto_prim(point_graph, alpha, axon=axon) 
        assert is_tree(pareto_tree2)
        mcost2, scost2 = graph_costs(pareto_tree2)
         
        mcosts1.append(mcost1)
        scosts1.append(scost1)
       
        mcosts2.append(mcost2)
        scosts2.append(scost2) 


    min_mcost = min(mcosts1 + mcosts2)
    max_mcost = max(mcosts1 + mcosts2)
    #min_scost = min(scosts1 + scosts2)
    min_scost = opt_scost
    max_scost = max(scosts1 + scosts2)

    normalize_scost = make_normalize_function(min_scost)
    normalize_mcost = make_normalize_function(min_mcost)

    '''
    norm_mcosts1 = pylab.array(mcosts1) / max_mcost
    norm_scosts1 = pylab.array(scosts1) / max_scost

    norm_mcosts2 = pylab.array(mcosts2) / max_mcost
    norm_scosts2 = pylab.array(scosts2) / max_scost
    '''

    norm_mcosts1 = normalize_mcost(pylab.array(mcosts1))
    norm_scosts1 = normalize_scost(pylab.array(scosts1))

    norm_mcosts2 = normalize_mcost(pylab.array(mcosts2))
    norm_scosts2 = normalize_scost(pylab.array(scosts2))

    pylab.figure()
    pylab.plot(norm_mcosts1, norm_scosts1, c = 'b')
    pylab.scatter(norm_mcosts1, norm_scosts1, c='b', label='steiner')
    
    pylab.plot(norm_mcosts2, norm_scosts2, c = 'k')
    pylab.scatter(norm_mcosts2, norm_scosts2, c='k', label='mst')
     
    pylab.xlabel('spanning tree cost')
    pylab.ylabel('satellite cost')
    
    #neural_mcost = normalize_mcost(mst_cost(G))
    #neural_scost  = normalize_scost(satellite_cost(G))
    neural_mcost, neural_scost = graph_costs(G)
    #neural_mcost /= max_mcost
    #neural_scost /= max_scost
    pylab.scatter([normalize_mcost(neural_mcost)], [normalize_scost(neural_scost)], c='r', marker='x', linewidths=15, label='neural')

    centroid_tree = centroid_mst(point_graph)
    #centroid_mcost = normalize_mcost(mst_cost(centroid_tree))
    #centroid_scost = normalize_scost(satellite_cost(centroid_tree))
    centroid_mcost, centroid_scost = graph_costs(centroid_tree)
    #centroid_mcost /= max_mcost
    #centroid_scost /= max_scost
    pylab.scatter([normalize_mcost(centroid_mcost)], [normalize_scost(centroid_scost)], c='g', marker='+', linewidths=15, label='centroid')

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
                            outdir = 'steiner_figs/%s/%s/%s/%s/%s' % (cell_type, species, region, lab, name)
                            outdir = outdir.replace(' ', '_')
                            os.system('mkdir -p %s' % outdir)
                            if len(os.listdir(outdir)) > 0:
                                continue
                            print species, lab, neuron
                            axon = i == 0
                            try:
                                pareto_plot_neuromorpho(G, name, cell_type,\
                                                        species, region,\
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
            outdir = 'steiner_imaris/' + subdir
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-min_nodes', type=int, default=MIN_NODES)
    parser.add_argument('-max_nodes', type=int, default=MAX_NODES)
    parser.add_argument('-a', '--algorithm', choices=['greedy', 'khuller'], default='greedy')
    parser.add_argument('-n', '--neuromorpho', action='store_true')
    parser.add_argument('-i', '--imaris', action='store_true')
    parser.add_argument('-b', '--neuron_builder', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')

    args = parser.parse_args()
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes
    algorithm = args.algorithm
    neuromorpho = args.neuromorpho
    imaris = args.imaris
    neuron_builder = args.neuron_builder
    debug = args.debug

    if debug:
        cell_type = 'principal_cell'
        species = 'giraffe'
        region = 'neocortex'
        lab = 'Jacobs'
        name = '185-2-24dk0'
        neuron = '185-2-24dk.CNG.swc'
        filename = 'datasets/%s/%s/%s/%s/%s' % (cell_type, species, region, lab, neuron)

        graphs = get_neuron_points(filename)
        G = graphs[0]

        pareto_plot_neuromorpho(G, name, cell_type, species, region, lab)

        return None

    if imaris:
        imaris_plots()
    if neuromorpho:
        neuromorpho_plots(min_nodes, max_nodes)
    if neuron_builder:
        neuron_builder_plots()

if __name__ == '__main__':
    main()
