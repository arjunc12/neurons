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
import math

SKIP_TYPES = ['Unknown_neurotransmitter', 'Not_reported']

VIZ_TREES = False
MIN_NODES = 0
MAX_NODES = 3000
NO_CONTINUE = False

LOG_PLOT = True
def ceil_power_of_10(n):
    exp = math.log(n, 10)
    exp = math.ceil(exp)
    return 10**exp

def floor_power_of_10(n):
    exp = math.log(n, 10)
    exp = math.ceil(exp)
    return 10**exp

def pareto_plot_neuromorpho(G, name, cell_type, species, region, lab,\
                            outdir='steiner_figs', output=True,\
                            viz_trees=VIZ_TREES, axon=False):
    assert G.number_of_nodes() > 0

    assert is_tree(G)
   
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
    opt_mcost = mst_cost(span_tree) / 2.0
 
# ---------------------------------------
    mcosts = []
    scosts = []
    costs = []

    delta = 0.01
    alphas = np.arange(delta, 1 + delta, delta)
    #alphas = [0.99]
    print "sorting neighbors"
    sort_neighbors(point_graph)

    khuller_comparisons = 0
    khuller_dominates = 0

    os.system("mkdir -p %s" % outdir)

# ---------------------------------------
    for i, alpha in enumerate(alphas):
        print "alpha", alpha
        pareto_tree = pareto_steiner(point_graph, alpha, axon=axon)
        #pareto_tree = pareto_steiner_sandbox(point_graph, alpha, axon=axon)
            
        assert is_tree(pareto_tree)
        mcost = mst_cost(pareto_tree)
        scost = satellite_cost(pareto_tree, relevant_nodes=point_graph.nodes())
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
       
        '''
        pareto_tree2 = pareto_khuller(point_graph, alpha, span_tree, sat_tree)
        assert is_tree(pareto_tree2)
        mcost2, scost2 = graph_costs(pareto_tree2)
        cost2 = pareto_cost(mcost=mcost2, scost=scost2, alpha=alpha)

        khuller_comparisons += 1
        if cost <= cost2:
            khuller_dominates += 1
        '''

# ---------------------------------------
    min_mcost = min(mcosts)
    max_mcost = max(mcosts)
    min_scost = min(scosts)
    max_scost = max(scosts)

    opt_mcost = min_mcost
    opt_scost = min_scost

    normalize_scost = make_normalize_function(opt_scost)
    normalize_mcost = make_normalize_function(opt_mcost)

    norm_mcosts = map(normalize_mcost, mcosts)
    norm_scosts = map(normalize_scost, scosts)

# ---------------------------------------
    
    neural_mcost, neural_scost = graph_costs(G)
    
    neural_dist, neural_index = pareto_dist(mcosts, scosts, neural_mcost,\
                                            neural_scost)
    neural_closem = mcosts[neural_index] 
    neural_closes = scosts[neural_index]
    neural_alpha = alphas[neural_index]
    
# ---------------------------------------
    centroid_tree = centroid_mst(point_graph)

    centroid_mcost, centroid_scost = graph_costs(centroid_tree)
    norm_centroid_mcost = normalize_mcost(centroid_mcost)
    norm_centroid_scost = normalize_scost(centroid_scost)
    
    centroid_dist, centroid_index = pareto_dist(mcosts, scosts, centroid_mcost, centroid_scost)
    centroid_closem = mcosts[centroid_index]
    centroid_closes = scosts[centroid_index]
    centroid_alpha = alphas[centroid_index]

# ---------------------------------------

    norm_neural_mcost = normalize_mcost(neural_mcost)
    norm_neural_scost = normalize_scost(neural_scost)


    norm_neural_dist, norm_neural_index = pareto_dist(norm_mcosts,\
                                                      norm_scosts,\
                                                      norm_neural_mcost,\
                                                      norm_neural_scost)
    norm_neural_alpha = alphas[norm_neural_index]

    norm_centroid_dist, norm_centroid_index = pareto_dist(norm_mcosts,\
                                                          norm_scosts,\
                                                          norm_centroid_mcost,\
                                                          norm_centroid_scost)
    norm_centroid_alpha = alphas[norm_centroid_index] 
        
# ---------------------------------------
    ntrials = 20
    successes = 0
    norm_successes = 0

    total_rand_dist = 0
    total_norm_rand_dist = 0
    rand_mcosts = []
    rand_scosts = []

    norm_rand_mcosts = []
    norm_rand_scosts = []
    for i in xrange(ntrials):
        rand_mst = random_mst(point_graph)
        
        rand_mcost, rand_scost = graph_costs(rand_mst)
        
        rand_mcosts.append(rand_mcost)
        rand_scosts.append(rand_scost)
        
        rand_dist, rand_index = pareto_dist(mcosts, scosts, rand_mcost, rand_scost)
        rand_closem = mcosts[rand_index]
        rand_closes = scosts[rand_index]
        rand_alpha = alphas[rand_index]
       
        norm_rand_mcost = normalize_mcost(rand_mcost)
        norm_rand_scost = normalize_scost(rand_scost)
        norm_rand_mcosts.append(norm_rand_mcost)
        norm_rand_scosts.append(norm_rand_scost)
        norm_rand_dist, norm_rand_index = pareto_dist(norm_mcosts,\
                                                      norm_scosts,\
                                                      norm_rand_mcost,\
                                                      norm_rand_scost)
        norm_rand_alpha = alphas[norm_rand_index]

        
        if rand_dist < neural_dist:
            successes += 1
        if norm_rand_dist < norm_neural_dist:
            norm_successes += 1

        total_rand_dist += rand_dist
        total_norm_rand_dist += norm_rand_dist

# ---------------------------------------
    pylab.figure()

    pylab.plot(mcosts, scosts, c = 'b')
    pylab.scatter(mcosts, scosts, c='b', label='greedy steiner')
   
    pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x',\
                   linewidths=15, label='neural') 
    
    pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+',
                   linewidths=15, label='centroid mst')
        
    pylab.scatter(rand_mcosts, rand_scosts, c='m', marker='o', label='random mst')
    
    #pylab.scatter(rand_mcosts, rand_scosts, c='m', marker='o', label='random mst')
     
    pylab.xlabel('steiner tree cost')
    pylab.ylabel('satellite cost')

    pylab.xlim(min_mcost - 10, max_mcost + 10)
    pylab.ylim(min_scost - 10, max_scost + 10)
    
    pylab.savefig('%s/pareto_front_%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

# ---------------------------------------
    if LOG_PLOT:
        pylab.figure()
        
        mcosts = pylab.log10(mcosts)
        scosts = pylab.log10(scosts)
        
        neural_mcost = pylab.log10(neural_mcost)
        neural_scost = pylab.log10(neural_scost)
        
        centroid_mcost = pylab.log10(centroid_mcost)
        centroid_scost = pylab.log10(centroid_scost)

        rand_mcosts = pylab.log10(rand_mcosts)
        rand_scosts = pylab.log10(rand_scosts)

        pylab.plot(mcosts, scosts, c = 'b')
        pylab.scatter(mcosts, scosts, c='b', label='greedy steiner')
       
        pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x',\
                       linewidths=15, label='neural') 
        
        pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+',
                       linewidths=15, label='centroid mst')
        
        pylab.scatter(rand_mcosts, rand_scosts, c='m', marker='o', label='random mst')
         
        pylab.xlabel('steiner tree cost')
        pylab.ylabel('satellite cost')
        
        pylab.savefig('%s/log-pareto_front_%s.pdf' % (outdir, name), format='pdf')
        pylab.close()
    
# ---------------------------------------
    pylab.figure()
    pylab.plot(norm_mcosts, norm_scosts, c = 'b')
    pylab.scatter(norm_mcosts, norm_scosts, c='b', label='greedy steiner')
    
    pylab.scatter([norm_neural_mcost], [norm_neural_scost], c='r', marker='x',\
                   linewidths=15, label='neural')
    
    pylab.scatter([norm_centroid_mcost], [norm_centroid_scost], c='g', marker='+',\
                   linewidths=15, label='centroid mst')
     
    pylab.scatter(norm_rand_mcosts, norm_rand_scosts, c='m', marker='o', label='random mst')
    
    pylab.xlabel('steiner tree cost')
    pylab.ylabel('satellite cost')
    
    pylab.savefig('%s/norm_pareto_front_%s.pdf' % (outdir, name), format='pdf')
    pylab.close()
    
# ---------------------------------------

    viz_tree(G, name + str('_neural'), outdir=outdir) 
    viz_tree(sat_tree, name + str('_sat'), outdir=outdir)
    viz_tree(span_tree, name + str('_mst'), outdir=outdir)
    viz_tree(centroid_tree, name + '_centroid', outdir=outdir) 

    mean_rand_dist = total_rand_dist / float(ntrials)
    mean_norm_rand_dist = total_norm_rand_dist / float(ntrials)
 
    def remove_spaces(string):
        return string.replace(' ', '')

    def remove_commas(string):
        return string.replace(',', '')
    f = open('pareto_steiner.csv', 'a')
    
    if output:
        write_items = [name, cell_type, species, region, lab, point_graph.number_of_nodes()]
        
        write_items.append(neural_alpha)
        write_items.append(norm_neural_alpha)
        
        write_items += [neural_dist, centroid_dist, mean_rand_dist]
        write_items += [norm_neural_dist, norm_centroid_dist, mean_norm_rand_dist]
        
        write_items.append(ntrials)
        write_items.append(successes)
        write_items.append(norm_successes)
        
        write_items = map(str, write_items)
        write_items = map(remove_commas, write_items)
        write_items = map(remove_spaces, write_items)
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

    mcosts3 = []
    scosts3 = []

    delta = 0.01
    alphas = np.arange(delta, 1, delta)
    print "sorting neighbors"
    sort_neighbors(point_graph)

    os.system("mkdir -p %s" % outdir)

    sat_tree = satellite_tree(point_graph)
    opt_scost = float(satellite_cost(sat_tree))
    span_tree = nx.minimum_spanning_tree(point_graph, weight='length')
    opt_mcost = float(mst_cost(span_tree))

    mst_comparisons = 0
    mst_dominates = 0

    khuller_comparisons = 0
    khuller_dominates = 0

    root = G.graph['root']

    for i, alpha in enumerate(alphas):
        print "alpha", alpha
        pareto_tree1 = pareto_steiner(point_graph, alpha, axon=axon)
        #pareto_tree1 = pareto_steiner_sandbox(point_graph, alpha, axon=axon)
        assert is_tree(pareto_tree1)
        mcost1 = float(mst_cost(pareto_tree1))
        scost1 = float(satellite_cost(pareto_tree1, relevant_nodes=point_graph.nodes()))
        cost1 = pareto_cost(mcost=mcost1, scost=scost1, alpha=alpha)
        
        if scost1 < opt_scost:
            print scost1, opt_scost, opt_scost - scost1
            print "dist error", dist_error(pareto_tree1)
            print "scost error", scost_error(pareto_tree1)
        assert scost1 >= opt_scost


        if False and (i % 5 == 0) and viz_trees:
            viz_tree(pareto_tree1, name + '-' + str(alpha), outdir=outdir)
        
        pareto_tree2 = pareto_prim(point_graph, alpha, axon=axon) 

        '''
        label_points(pareto_tree1)
        label_points(pareto_tree2)
        viz_tree(pareto_tree1, name + '_steiner' + str(alpha), outdir='debug_steiner')
        viz_tree(pareto_tree2, name + '_mst' + str(alpha), outdir='debug_steiner')
        '''

        assert is_tree(pareto_tree2)
        mcost2, scost2 = graph_costs(pareto_tree2)
        cost2 = pareto_cost(mcost=mcost2, scost=scost2, alpha=alpha)
        
        pareto_tree3 = pareto_khuller(point_graph, alpha, span_tree, sat_tree)
        assert is_tree(pareto_tree3)
        mcost3, scost3 = graph_costs(pareto_tree3)
        cost3 = pareto_cost(mcost=mcost3, scost=scost3, alpha=alpha)

        assert pareto_tree2.number_of_nodes() == pareto_tree3.number_of_nodes()

        '''
        print "steiner", mcost1, scost1, cost1
        print "prim   ", mcost2, scost2, cost2
        print "khuller", mcost3, scost3, cost3
        '''

        mcosts1.append(mcost1)
        scosts1.append(scost1)
       
        mcosts2.append(mcost2)
        scosts2.append(scost2) 

        mcosts3.append(mcost3)
        scosts3.append(scost3)
 
        mst_comparisons += 1
        if cost1 <= cost2:
            mst_dominates += 1

        khuller_comparisons += 1
        if cost1 <= cost3:
            khuller_dominates += 1

    outfile = open('steiner_imaris.csv', 'a')
    write_items = [mst_comparisons, mst_dominates,\
                   khuller_comparisons, khuller_dominates]
    write_items = map(str, write_items)
    write_items = ', '.join(write_items)
    outfile.write('%s\n' % write_items)

    outfile.write('khuller dominated by steiner, ' + str(prop_dominated(mcosts1, scosts1, mcosts3, scosts3)) + '\n')
    outfile.write('steiner dominated by khuller, ' + str(prop_dominated(mcosts3, scosts3, mcosts1, scosts1)) + '\n')
    outfile.write('mst dominated by steiner, ' + str(prop_dominated(mcosts1, scosts1, mcosts2, scosts2)) + '\n')
    outfile.write('steiner dominated by mst, ' + str(prop_dominated(mcosts2, scosts2, mcosts1, scosts1)) + '\n')

    outfile.close()
    
    min_mcost = min(mcosts1 + mcosts2 + mcosts3)
    max_mcost = max(mcosts1 + mcosts2 + mcosts3)
    #min_scost = min(scosts1 + scosts2)
    min_scost = opt_scost
    max_scost = max(scosts1 + scosts2 + scosts3)

    normalize_scost = make_normalize_function(min_scost)
    normalize_mcost = make_normalize_function(min_mcost)


    norm_mcosts1 = normalize_mcost(pylab.array(mcosts1))
    norm_scosts1 = normalize_scost(pylab.array(scosts1))

    norm_mcosts2 = normalize_mcost(pylab.array(mcosts2))
    norm_scosts2 = normalize_scost(pylab.array(scosts2))

    norm_mcosts3 = normalize_mcost(pylab.array(mcosts3))
    norm_scosts3 = normalize_scost(pylab.array(scosts3))

    '''
    pylab.figure()
    pylab.plot(norm_mcosts1, norm_scosts1, c = 'r')
    pylab.scatter(norm_mcosts1, norm_scosts1, c='r', label='steiner')
    
    pylab.plot(norm_mcosts2, norm_scosts2, c = 'b')
    pylab.scatter(norm_mcosts2, norm_scosts2, c='b', label='mst')

    pylab.plot(norm_mcosts3, norm_scosts3, c='k')
    pylab.scatter(norm_mcosts3, norm_scosts3, c='k', label='khuller')
     
    pylab.xlabel('spanning tree cost')
    pylab.ylabel('satellite cost')
    '''

    pylab.figure()
    pylab.plot(mcosts1, scosts1, c = 'b')
    pylab.scatter(mcosts1, scosts1, c='b', label='steiner')
    
    #pylab.plot(mcosts2, scosts2, c = 'r')
    #pylab.scatter(mcosts2, scosts2, c='r', label='mst')

    #pylab.plot(mcosts3, scosts3, c='k')
    #pylab.scatter(mcosts3, scosts3, c='k', label='khuller')
     
    pylab.xlabel('spanning tree cost')
    pylab.ylabel('satellite cost')

    neural_mcost, neural_scost = graph_costs(G)
    
    norm_neural_mcost = normalize_mcost(neural_mcost)
    norm_neural_Scost = normalize_scost(neural_scost)
    
    #pylab.scatter([norm_neural_mcost], [norm_neural_scost], c='r', marker='x', linewidths=15, label='neural')
    pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x', linewidths=15, label='neural')

    centroid_tree = centroid_mst(point_graph)
    centroid_mcost, centroid_scost = graph_costs(centroid_tree)

    norm_centroid_mcost = normalize_mcost(centroid_mcost)
    norm_centroid_scost = normalize_scost(centroid_scost)

    #pylab.scatter([norm_centroid_mcost)], [norm_centroid_scost], c='g', marker='+', linewidths=15, label='centroid')
    #pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+', linewidths=15, label='centroid')

    pylab.savefig('%s/pareto_front_%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

    viz_tree(G, name + str('_neural'), outdir=outdir) 
    viz_tree(sat_tree, name + str('_sat'), outdir=outdir)
    viz_tree(span_tree, name + str('_mst'), outdir=outdir)
    viz_tree(centroid_tree, name + '_centroid', outdir=outdir) 

def neuromorpho_plots(min_nodes=MIN_NODES, max_nodes=MAX_NODES, cell_types=None,\
                      animal_species=None, regions=None):
    #directory = 'neuromorpho'
    directory = 'datasets'
    i = 0
    for cell_type in os.listdir(directory):
        if cell_type in SKIP_TYPES:
            continue
        if cell_types != None and cell_type not in cell_types:
            continue
        for species in os.listdir(directory + '/' + cell_type):
            if animal_species != None and species not in animal_species:
                continue
            for region in os.listdir(directory + '/' + cell_type + '/' + species):
                if regions != None and region not in regions:
                    continue
                for lab in os.listdir(directory + "/" + cell_type + '/' + species+ '/' + region):
                    for neuron in os.listdir(directory + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = directory + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron
                        
                        if neuron[-8:] != ".CNG.swc": 
                            continue
                        
                        try:
                            graphs = get_neuron_points(filename)
                        except AssertionError as e:
                            continue

                        for i, G in enumerate(graphs):
                            if G == None:
                                continue
                            if not (min_nodes <= G.number_of_nodes() <= max_nodes):
                                print "wrong nodes", G.number_of_nodes()
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
    os.system('rm -f steiner_imaris.csv')
    i = 0
    for subdir in os.listdir('imaris'):
        if os.path.isdir('imaris/' + subdir):
            i += 1
            if i != 1:
                pass
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
    #parser.add_argument('-a', '--algorithm', choices=['greedy', 'khuller'], default='greedy')
    parser.add_argument('-n', '--neuromorpho', action='store_true')
    parser.add_argument('-i', '--imaris', action='store_true')
    parser.add_argument('-b', '--neuron_builder', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-c', '--cell_types', nargs='+', default=None)
    parser.add_argument('-s', '-a', '--species', '--animal_species,', nargs='+',\
                        default=None, dest='animal_species')
    parser.add_argument('-r', '--regions', nargs='+', default=None)

    args = parser.parse_args()
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes
    #algorithm = args.algorithm
    neuromorpho = args.neuromorpho
    imaris = args.imaris
    neuron_builder = args.neuron_builder
    debug = args.debug
    cell_types = args.cell_types
    animal_species = args.animal_species
    regions = args.regions

    if debug:
        cell_type = 'bipolar'
        species = 'human'
        region = 'retina'
        lab = 'kantor'
        name = 'humret_FMB_40x_5'
        neuron = name + '.CNG.swc'
        graph_number = 1
        name += str(graph_number)
        filename = 'datasets/%s/%s/%s/%s/%s' % (cell_type, species, region, lab, neuron)

        graphs = get_neuron_points(filename)
        G = graphs[graph_number]

        pareto_plot_neuromorpho(G, name, cell_type, species, region, lab)

        return None

    if imaris:
        imaris_plots()
    if neuromorpho:
        neuromorpho_plots(min_nodes, max_nodes, cell_types, animal_species, regions)
    if neuron_builder:
        neuron_builder_plots()

if __name__ == '__main__':
    main()
