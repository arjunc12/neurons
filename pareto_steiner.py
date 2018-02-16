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
from random_graphs import random_mst, barabasi_tree
import math
import seaborn as sns
from collections import defaultdict
from numpy.random import permutation
import pandas as pd

SKIP_TYPES = ['Unknown_neurotransmitter', 'Not_reported']

VIZ_TREES = True
VIZ_TREES_NEUROMORPHO = False
VIZ_TREES_IMARIS = True
MIN_NODES = 0
MAX_NODES = 3000
NO_CONTINUE = False

LOG_PLOT = True

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'
PARETO_FIGS_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_figs'
PARETO_FRONT_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts'
PARETO_OUTPUT_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_steiner_temp'

# DIST_FUNC = pareto_dist_l2
DIST_FUNC = pareto_dist_scale

NEUROMORPHO_PADDING = 500

NEURON_TYPES = {0 : 'axon', 1 : 'basal dendrite', 2: 'apical dendrite',\
                3: 'truncated axon'}

def ceil_power_of_10(n):
    exp = math.log(n, 10)
    exp = math.ceil(exp)
    return 10**exp

def floor_power_of_10(n):
    exp = math.log(n, 10)
    exp = math.ceil(exp)
    return 10**exp

def pareto_plot_neuromorpho(G, neuron_name, neuron_type,\
                            front_dir=PARETO_FRONT_DIR,\
                            figs_dir=PARETO_FRONT_DIR,\
                            viz_trees=VIZ_TREES_NEUROMORPHO):
    assert G.number_of_nodes() > 0
    assert is_tree(G)
   
    #print "making graph"
    synapses = G.graph['synapses']
    points = synapses + [G.graph['root']]
    point_graph = G.subgraph(points)
    point_graph = complete_graph(point_graph)

    sat_tree = satellite_tree(point_graph)
    sat_scost = satellite_cost(sat_tree)
    sat_mcost = mst_cost(sat_tree)
    span_tree = nx.minimum_spanning_tree(point_graph, weight='length') 

# ---------------------------------------
    pareto_front_fname = '%s/pareto_front.csv' % front_dir
    first_time = not os.path.exists(pareto_front_fname)

    tree_costs_fname = '%s/tree_costs.csv' % front_dir
   
    models_fname = '%s/models_%s.csv' % (output_dir, neuron_name)
    
    output_fname = '%s/pareto_steiner_%s.csv' %  (output_dir, neuron_name)

# ---------------------------------------
    alphas = None
    mcosts = None
    scosts = None

# ---------------------------------------
    if (not first_time) and (not viz_trees):
        df = pd.read_csv('%s/pareto_front.csv' % front_dir, skipinitialspace=True)
        alphas = pylab.array(df['alpha'])
        mcosts = pylab.array(df['mcost'])
        scosts = pylab.array(df['scost'])
    else:
        delta = 0.01
        alphas = np.arange(delta, 1 + delta, delta)
        mcosts = []
        scosts = []
        costs = []

        print "sorting neighbors"
        sort_neighbors(point_graph)

        front_lines = []
        front_lines.append('alpha, mcost, scost\n')
        front_lines.append('%f, %f, %f\n' % (0, sat_mcost, sat_scost))
        axon = neuron_type == 'axon'
        for i, alpha in enumerate(alphas):
            print alpha
            pareto_tree = pareto_steiner(point_graph, alpha, axon=axon)
            #pareto_tree = pareto_steiner_sandbox(point_graph, alpha, axon=axon)
                
            #assert is_tree(pareto_tree)
            mcost = mst_cost(pareto_tree)
            scost = satellite_cost(pareto_tree, relevant_nodes=point_graph.nodes())
            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
            '''
            if scost < sat_scost:
                print scost, opt_scost, opt_scost - scost
                print "dist error", dist_error(pareto_tree)
                print "scost error", scost_error(pareto_tree)
            assert scost >= sat_scost
            check_dists(pareto_tree)
            '''

            if (i % 5 == 4) and viz_trees:
                viz_tree(pareto_tree, '%s-%0.2f' % ('alpha', alpha), outdir=figs_dir)
            
            mcosts.append(mcost)
            scosts.append(scost)

            front_lines.append('%f, %f, %f\n' % (alpha, mcost, scost))

        front_fname = '%s/pareto_front.csv' % front_dir
        with open(front_fname, 'w') as front_file:
            front_file.writelines(front_lines)
    pass

def pareto_front_neuromorpho(G, neuron_name, neuron_type,\
                             front_dir=PARETO_FRONT_DIR,\
                             output_dir=PARETO_OUTPUT_DIR,\
                             figs_dir=PARETO_FIGS_DIR, output=True,\
                             viz_trees=VIZ_TREES_NEUROMORPHO):
    
    assert G.number_of_nodes() > 0
    assert is_tree(G)
   
    #print "making graph"
    synapses = G.graph['synapses']
    points = synapses + [G.graph['root']]
    point_graph = G.subgraph(points)
    point_graph = complete_graph(point_graph)

    sat_tree = satellite_tree(point_graph)
    sat_scost = satellite_cost(sat_tree)
    sat_mcost = mst_cost(sat_tree)
    span_tree = nx.minimum_spanning_tree(point_graph, weight='length') 

# ---------------------------------------
    pareto_front_fname = '%s/pareto_front.csv' % front_dir
    first_time = not os.path.exists(pareto_front_fname)

    tree_costs_fname = '%s/tree_costs.csv' % front_dir
   
    models_fname = '%s/models_%s.csv' % (output_dir, neuron_name)
    
    output_fname = '%s/pareto_steiner_%s.csv' %  (output_dir, neuron_name)

# ---------------------------------------
    alphas = None
    mcosts = None
    scosts = None

# ---------------------------------------
    if (not first_time) and (not viz_trees):
        df = pd.read_csv('%s/pareto_front.csv' % front_dir, skipinitialspace=True)
        alphas = pylab.array(df['alpha'])
        mcosts = pylab.array(df['mcost'])
        scosts = pylab.array(df['scost'])
    else:
        delta = 0.01
        alphas = np.arange(delta, 1 + delta, delta)
        mcosts = []
        scosts = []
        costs = []

        print "sorting neighbors"
        sort_neighbors(point_graph)

        front_lines = []
        front_lines.append('alpha, mcost, scost\n')
        front_lines.append('%f, %f, %f\n' % (0, sat_mcost, sat_scost))
        axon = neuron_type == 'axon'
        for i, alpha in enumerate(alphas):
            print alpha
            pareto_tree = pareto_steiner(point_graph, alpha, axon=axon)
            #pareto_tree = pareto_steiner_sandbox(point_graph, alpha, axon=axon)
                
            #assert is_tree(pareto_tree)
            mcost = mst_cost(pareto_tree)
            scost = satellite_cost(pareto_tree, relevant_nodes=point_graph.nodes())
            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
            '''
            if scost < sat_scost:
                print scost, opt_scost, opt_scost - scost
                print "dist error", dist_error(pareto_tree)
                print "scost error", scost_error(pareto_tree)
            assert scost >= sat_scost
            check_dists(pareto_tree)
            '''

            if (i % 5 == 4) and viz_trees:
                viz_tree(pareto_tree, '%s-%0.2f' % ('alpha', alpha), outdir=figs_dir)
            
            mcosts.append(mcost)
            scosts.append(scost)

            front_lines.append('%f, %f, %f\n' % (alpha, mcost, scost))

        front_fname = '%s/pareto_front.csv' % front_dir
        with open(front_fname, 'w') as front_file:
            front_file.writelines(front_lines)
   
# ---------------------------------------
    
    #neural_mcost, neural_scost = graph_costs(G)
    neural_mcost = mst_cost(G)
    neural_scost = satellite_cost(G, relevant_nodes=point_graph.nodes())
    
    neural_dist, neural_index = DIST_FUNC(mcosts, scosts, neural_mcost,\
                                          neural_scost)
    neural_closem = mcosts[neural_index] 
    neural_closes = scosts[neural_index]
    neural_alpha = alphas[neural_index]

    
# ---------------------------------------
    centroid_tree = centroid_mst(point_graph)

    centroid_mcost, centroid_scost = graph_costs(centroid_tree)
    
    centroid_dist, centroid_index = DIST_FUNC(mcosts, scosts,\
                                              centroid_mcost,\
                                              centroid_scost)
    centroid_closem = mcosts[centroid_index]
    centroid_closes = scosts[centroid_index]
    centroid_alpha = alphas[centroid_index]


# ---------------------------------------
    centroid_success = int(centroid_dist <= neural_dist)
    centroid_ratio = centroid_dist / neural_dist
    if first_time:
        with open(models_fname, 'w') as models_file:
            models_file.write('%s, %s, %s, %f,,\n' % (neuron_name,\
                                                      neuron_type,\
                                                      'neural',\
                                                      neural_dist))
            models_file.write('%s, %s, %s, %f, %d, %f\n' % (neuron_name,\
                                                            neuron_type,\
                                                            'centroid',\
                                                            centroid_dist,\
                                                            centroid_success,\
                                                            centroid_ratio))

# ---------------------------------------
    if first_time:
        with open(tree_costs_fname, 'w') as tree_costs_file:
            tree_costs_file.write('tree, mcost, scost\n')
            tree_costs_file.write('%s, %f, %f\n' % ('neural',neural_mcost,\
                                                    neural_scost))
            tree_costs_file.write('%s, %f, %f\n' % ('centroid', centroid_mcost,\
                                                    centroid_scost))
# ---------------------------------------
    random_trials = 20

    for i in xrange(random_trials):
        rand_mst = random_mst(point_graph)
        rand_mcost, rand_scost = graph_costs(rand_mst)
        rand_dist, rand_index = DIST_FUNC(mcosts, scosts, rand_mcost,\
                                          rand_scost)
        rand_success = int(rand_dist <= neural_dist)
        rand_ratio = rand_dist / neural_dist

        barabasi_mst = barabasi_tree(point_graph)
        barabasi_mcost, barabasi_scost = graph_costs(barabasi_mst)
        barabasi_dist, barabasi_index = DIST_FUNC(mcosts, scosts,\
                                                  barabasi_mcost,\
                                                  barabasi_scost)
        barabasi_success = int(barabasi_dist <= neural_dist)
        barabasi_ratio = barabasi_dist / neural_dist

        with open(tree_costs_fname, 'a') as tree_costs_file:
            tree_costs_file.write('%s, %f, %f\n' % ('random', rand_mcost,\
                                                     rand_scost))
            tree_costs_file.write('%s, %f, %f\n' % ('barabasi',\
                                                     barabasi_mcost,\
                                                     barabasi_scost))

        with open(models_fname, 'a') as models_file:
            models_file.write('%s, %s, %s, %f, %d, %f\n' % (neuron_name,\
                                                            neuron_type,\
                                                            'random',\
                                                            rand_dist,\
                                                            rand_success,\
                                                            rand_ratio))
            
            models_file.write('%s, %s, %s, %f, %d, %f\n' % (neuron_name,\
                                                            neuron_type,\
                                                            'barabasi',\
                                                            barabasi_dist,\
                                                            barabasi_success,\
                                                            barabasi_ratio))

# ---------------------------------------
    def remove_spaces(string):
        return string.replace(' ', '')

    def remove_commas(string):
        return string.replace(',', '')
     
    if output and first_time:
        write_items = [neuron_name, neuron_type, point_graph.number_of_nodes()]
        
        write_items.append(neural_alpha)
        
        write_items = map(str, write_items)
        write_items = map(remove_commas, write_items)
        #write_items = map(remove_spaces, write_items)
        write_items = ', '.join(write_items)
        with open(output_fname, 'w') as output_file:
            output_file.write('%s\n' % write_items)

# ---------------------------------------
    pylab.figure()
    sns.set()

    pylab.plot(mcosts, scosts, c = 'b')
    pylab.scatter(mcosts, scosts, c='b', label='greedy steiner')
   
    pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x',\
                   linewidths=60, label='neural', s=175) 
     
    pylab.xlabel('wiring cost')
    pylab.ylabel('conduction delay')
    
    pylab.savefig('%s/pareto_front.pdf' % figs_dir, format='pdf')
    pylab.close()

    
# ---------------------------------------

    for u in G.nodes():
        for H in [G, sat_tree, span_tree, centroid_tree]:
            assert H.has_node(u)
            if u == G.graph['root']:
                H.node[u]['label'] = 'root'
            else:
                H.node[u]['label'] = 'synapse'
    if viz_trees:
        viz_tree(G, 'neural', outdir=figs_dir) 
        viz_tree(sat_tree, 'sat', outdir=figs_dir)
        viz_tree(span_tree, 'mst', outdir=figs_dir)
        viz_tree(centroid_tree, 'centroid', outdir=figs_dir) 
 
def pareto_plot_imaris(G, name, outdir='steiner_imaris', viz_trees=VIZ_TREES_IMARIS, axon=False):
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
    alpas = np.append(alphas, [1])
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

    costs_file = open('steiner_imaris.csv', 'a')

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

        costs_file.write('%s, %f, %f, %f\n' % (name, alpha, mcost1, scost1))

        if (i % 5 == 4) and viz_trees:
            viz_tree(pareto_tree1, '%s-%0.2f' % (name, alpha), outdir=outdir)
            #viz_tree(pareto_tree1, name + '-' + str(alpha), outdir=outdir)
        
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

    costs_file.close()

    '''
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
    '''
    
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
    sns.set()
    pylab.plot(mcosts1, scosts1, c = 'b')
    pylab.scatter(mcosts1, scosts1, c='b', label='steiner')
    
    #pylab.plot(mcosts2, scosts2, c = 'r')
    #pylab.scatter(mcosts2, scosts2, c='r', label='mst')

    #pylab.plot(mcosts3, scosts3, c='k')
    #pylab.scatter(mcosts3, scosts3, c='k', label='khuller')
     
    #pylab.xlabel('steiner tree cost')
    #pylab.ylabel('satellite cost')
    pylab.xlabel('wiring cost')
    pylab.ylabel('conduction delay')
    

    neural_mcost, neural_scost = graph_costs(G)
    
    norm_neural_mcost = normalize_mcost(neural_mcost)
    norm_neural_Scost = normalize_scost(neural_scost)
    
    #pylab.scatter([norm_neural_mcost], [norm_neural_scost], c='r', marker='x', linewidths=15, label='neural')
    pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x',\
                  linewidths=60, label='neural', s=175)
    
    pylab.legend()

    pylab.savefig('%s/pareto_front_%s.pdf' % (outdir, name), format='pdf')
    #pylab.close()

    centroid_tree = centroid_mst(point_graph)
    centroid_mcost, centroid_scost = graph_costs(centroid_tree)

    norm_centroid_mcost = normalize_mcost(centroid_mcost)
    norm_centroid_scost = normalize_scost(centroid_scost)

    #pylab.scatter([norm_centroid_mcost)], [norm_centroid_scost], c='g', marker='+', linewidths=15, label='centroid')
    pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+',\
                  linewidths=60, label='centroid', s=175)

    neural_dist, neural_index = DIST_FUNC(mcosts1, scosts1, neural_mcost,\
                                          neural_scost)
    
    centroid_dist, centroid_index = DIST_FUNC(mcosts1, scosts1, centroid_mcost,\
                                              centroid_scost)
    
    for scale_factor, color in zip([neural_dist, centroid_dist], ['r', 'g']):
        x = scale_factor * pylab.array(mcosts1)
        y = scale_factor * pylab.array(scosts1)
        pylab.plot(x, y, c=color, linestyle='-')
        pylab.scatter(x, y, c=color, label='s = %0.2f' % scale_factor)

    pylab.legend()

    pylab.savefig('%s/pareto_front_scaled_%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

    viz_tree(G, name + str('_neural'), outdir=outdir) 
    viz_tree(sat_tree, name + str('_sat'), outdir=outdir)
    viz_tree(span_tree, name + str('_mst'), outdir=outdir)
    viz_tree(centroid_tree, name + '_centroid', outdir=outdir) 

def neuromorpho_plots(min_nodes=MIN_NODES, max_nodes=MAX_NODES, cell_types=None,\
                      animal_species=None, regions=None, labs=None, names=None):
    #directory = 'neuromorpho'

    for cell_type in permutation(os.listdir(DATASETS_DIR)):
        if cell_type in SKIP_TYPES:
            continue
        if cell_types != None and cell_type not in cell_types:
            continue
        for species in os.listdir(DATASETS_DIR + '/' + cell_type):
            if animal_species != None and species not in animal_species:
                continue
            for region in os.listdir(DATASETS_DIR + '/' + cell_type + '/' + species):
                if regions != None and region not in regions:
                    continue
                for lab in os.listdir(DATASETS_DIR + "/" + cell_type + '/' + species+ '/' + region):
                    if labs != None and lab not in labs:
                        continue
                    for neuron_file in os.listdir(DATASETS_DIR + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = DATASETS_DIR + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron_file
                        
                        neuron_name = neuron_file[:-8]
                        if names != None and neuron_name not in names:
                            continue

                        if neuron_file[-8:] != ".CNG.swc": 
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
                            
                            neuron_type = NEURON_TYPES[i]

                            front_dir = '%s/%s/%s' % (PARETO_FRONT_DIR, neuron_name, neuron_type)
                            front_dir = front_dir.replace(' ', '_')
                            os.system('mkdir -p %s' % front_dir)

                            output_dir = PARETO_OUTPUT_DIR
                            output_dir = output_dir.replace(' ', '_')
                            os.system('mkdir -p %s' % output_dir)

                            figs_dir = '%s/%s/%s/%s/%s/%s/%s' % (PARETO_FIGS_DIR,\
                                                                 cell_type,\
                                                                 species,\
                                                                 region, lab,\
                                                                 neuron_name,\
                                                                 neuron_type)

                            figs_dir = figs_dir.replace(' ', '_')
                            os.system('mkdir -p %s' % figs_dir)

                            try:
                                #H = add_synapses(G)
                                H = G.copy()
                                H.graph['synapses'] = []
                                for u in H.nodes_iter():
                                    if u != H.graph['root']:
                                        H.graph['synapses'].append(u)
                                
                                pareto_front_neuromorpho(H, neuron_name, neuron_type,\
                                                         front_dir=front_dir,\
                                                         output_dir=output_dir,\
                                                         figs_dir=figs_dir)
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

            neural_dist, neural_index = DIST_FUNC(mcosts, scosts, neural_mcost, neural_scost)
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

def boutons_plots():
    groups = defaultdict(list)
    for fname in os.listdir('boutons/swc_files'):
        group = fname[:-5]
        num = fname[-5]
        groups[group].append(num)

    all_edge_lengths = []
    for group in groups:
        G = None
        for num in groups[group]:
            fname = 'boutons/swc_files/' + group + num + '.swc'
            graphs = get_neuron_points(fname)
            H = graphs[0]
            if G == None:
                G = H
            else:
                G = nx.disjoint_union(G, H)
        viz_tree(G, group, 'boutons/drawings')
        edge_lengths = []
        for u, v in G.edges_iter():
            edge_lengths.append(G[u][v]['length'])
        print "group", group
        print pylab.mean(edge_lengths)
        pylab.figure()
        pylab.hist(edge_lengths)
        pylab.savefig('boutons/histograms/edge_lengths_%s.pdf' % group, format='pdf')
        pylab.close()
        all_edge_lengths += edge_lengths

    print "grand average"
    print pylab.mean(all_edge_lengths)
    pylab.figure()
    pylab.hist(all_edge_lengths)
    pylab.savefig('boutons/histograms/edge_lengths_all.pdf', format='pdf')
    pylab.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-min_nodes', type=int, default=MIN_NODES)
    parser.add_argument('-max_nodes', type=int, default=MAX_NODES)
    #parser.add_argument('-a', '--algorithm', choices=['greedy', 'khuller'], default='greedy')
    parser.add_argument('-n', '--neuromorpho', action='store_true')
    parer.add_argument('-np', '--neuromorpho_plot', action='store_true')
    parser.add_argument('-v', '--viz_trees', action='store_true')
    parser.add_argument('-i', '--imaris', action='store_true')
    parser.add_argument('-b', '--neuron_builder', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-c', '--cell_types', nargs='+', default=None)
    parser.add_argument('-s', '-a', '--species', '--animal_species,', nargs='+',\
                        default=None, dest='animal_species')
    parser.add_argument('-l', '--labs', nargs='+', default=None, dest='labs')
    parser.add_argument('-na', '--names', nargs='+', default=None, dest='names')
    parser.add_argument('-r', '--regions', nargs='+', default=None)
    parser.add_argument('-bt', '--boutons', action='store_true')

    args = parser.parse_args()
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes
    #algorithm = args.algorithm
    neuromorpho = args.neuromorpho
    neuromorpho_plot = args.neuromorpho_plot
    viz_trees = args.viz_trees
    imaris = args.imaris
    neuron_builder = args.neuron_builder
    debug = args.debug
    cell_types = args.cell_types
    animal_species = args.animal_species
    regions = args.regions
    labs = args.labs
    names = args.names
    boutons = args.boutons

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
        neuromorpho_plots(min_nodes, max_nodes, cell_types, animal_species, regions, labs, names)
    if neuron_builder:
        neuron_builder_plots()
    if boutons:
        boutons_plots()

if __name__ == '__main__':
    main()
