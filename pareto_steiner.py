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
from neuron_builder import build_neuron_snider
from random_graphs import random_mst, barabasi_tree
import math
import seaborn as sns
from collections import defaultdict
from numpy.random import permutation
import pandas as pd

SKIP_TYPES = ['Unknown_neurotransmitter', 'Not_reported']

VIZ_TREES = True
VIZ_TREES_NEUROMORPHO = False
VIZ_TREES_IMARIS = False
MIN_NODES = 0
MAX_NODES = 3000
NO_CONTINUE = False

LOG_PLOT = True

DATA_DRIVE = '/iblsn/data/Arjun/neurons'

NEUROMORPHO_DATASETS_DIR = '%s/datasets' % DATA_DRIVE
NEUROMORPHO_FIGS_DIR = '%s/pareto_steiner_output/steiner_figs' % DATA_DRIVE
NEUROMORPHO_OUTPUT_DIR = '%s/pareto_steiner_output' % DATA_DRIVE
NEUROMORPHO_FRONTS_DIR = '%s/pareto_fronts' % NEUROMORPHO_OUTPUT_DIR
NEUROMORPHO_TEMP_DIR = '%s/pareto_steiner_temp' % NEUROMORPHO_OUTPUT_DIR

IMARIS_DRIVE = '%s/imaris' % DATA_DRIVE
IMARIS_TRACINGS_DIR = '%s/tracings' % IMARIS_DRIVE
IMARIS_FIGS_DIR = '%s/figs' % IMARIS_DRIVE
IMARIS_FRONTS_DIR = '%s/fronts' % IMARIS_DRIVE

# DIST_FUNC = pareto_dist_l2
DIST_FUNC = pareto_dist_scale

NEUROMORPHO_PADDING = 500

NEURON_TYPES = {0 : 'axon', 1 : 'basal dendrite', 2: 'apical dendrite',\
                3: 'truncated axon'}

COLORS = {'neural' : 'r', 'centroid' : 'g', 'random' : 'm', 'barabasi' : 'c'}

PLOT_TREES = ['neural', 'centroid']

def ceil_power_of_10(n):
    exp = math.log(n, 10)
    exp = math.ceil(exp)
    return 10**exp

def floor_power_of_10(n):
    exp = math.log(n, 10)
    exp = math.ceil(exp)
    return 10**exp

def pareto_plot(fronts_dir, figs_dir):
    pareto_front = pd.read_csv('%s/pareto_front.csv' % fronts_dir,\
                               skipinitialspace=True)
    mcosts = pareto_front['mcost']
    scosts = pareto_front['scost']

    tree_costs = pd.read_csv('%s/tree_costs.csv' % fronts_dir,\
                             skipinitialspace=True)
    
    pylab.figure()
    sns.set()

    pylab.plot(mcosts, scosts, c = 'b', label='_nolegend_')
    pylab.scatter(mcosts, scosts, c='b', label='pareto front')
  
    for tree, costs in tree_costs.groupby('tree'):
        if tree in PLOT_TREES:
            pylab.scatter(costs['mcost'], costs['scost'], label=tree,\
                          marker='x', s=175, c=COLORS[tree])
     
    pylab.xlabel('wiring cost')
    pylab.ylabel('conduction delay')
    pylab.legend()
    
    pylab.savefig('%s/pareto_front.pdf' % figs_dir, format='pdf')
    pylab.close()

def pareto_tree_costs(G, point_graph, axon=False, viz_trees=False, figs_dir=None,\
                      sandbox=False):
    delta = 0.01
    alphas = np.arange(delta, 1 + delta, delta)
    mcosts = []
    scosts = []

    print "sorting neighbors"
    sort_neighbors(point_graph)

    pareto_func = pareto_steiner
    if sandbox:
        pareto_func = pareto_steiner_sandbox
    for i, alpha in enumerate(alphas):
        print alpha
        #pareto_tree = pareto_steiner(point_graph, alpha, axon=axon)
        pareto_tree = pareto_func(point_graph, alpha, axon=axon)
        mcost = mst_cost(pareto_tree)
        scost = satellite_cost(pareto_tree, relevant_nodes=point_graph.nodes())
        cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)

        if (i % 5 == 4) and viz_trees:
            assert figs_dir != None
            viz_tree(pareto_tree, '%s-%0.2f' % ('alpha', alpha), outdir=figs_dir)
        
        mcosts.append(mcost)
        scosts.append(scost)

    return alphas, mcosts, scosts

def pareto_front(G, point_graph, neuron_name, neuron_type,\
                 fronts_dir=NEUROMORPHO_FRONTS_DIR, figs_dir=NEUROMORPHO_FIGS_DIR,\
                 viz_trees=VIZ_TREES_NEUROMORPHO, sandbox=False):

    sat_tree = satellite_tree(point_graph)
    sat_scost = satellite_cost(sat_tree)
    sat_mcost = mst_cost(sat_tree)
    span_tree = nx.minimum_spanning_tree(point_graph, weight='length') 

# ---------------------------------------
    pareto_front_fname = '%s/pareto_front.csv' % fronts_dir
    first_time = not os.path.exists(pareto_front_fname)

# ---------------------------------------
    alphas = None
    mcosts = None
    scosts = None

# ---------------------------------------
    if (not first_time) and (not viz_trees):
        df = pd.read_csv('%s/pareto_front.csv' % fronts_dir, skipinitialspace=True)
        alphas = pylab.array(df['alpha'])
        mcosts = pylab.array(df['mcost'])
        scosts = pylab.array(df['scost'])
    else:
        front_lines = []
        front_lines.append('alpha, mcost, scost\n')
        front_lines.append('%f, %f, %f\n' % (0, sat_mcost, sat_scost))
        axon = neuron_type == 'axon'

        alphas, mcosts, scosts = pareto_tree_costs(G, point_graph, axon,\
                                                   viz_trees, figs_dir=figs_dir,\
                                                   sandbox=sandbox)

        for i in xrange(len(alphas)):
            alpha = alphas[i]
            mcost = mcosts[i]
            scost = scosts[i]
            front_lines.append('%f, %f, %f\n' % (alpha, mcost, scost))

        front_fname = '%s/pareto_front.csv' % fronts_dir
        with open(front_fname, 'w') as front_file:
            front_file.writelines(front_lines)
    
# ---------------------------------------
    for u in point_graph.nodes():
        for H in [G, sat_tree, span_tree]:
            assert H.has_node(u)
            if u == G.graph['root']:
                H.node[u]['label'] = 'root'
            else:
                H.node[u]['label'] = 'synapse'
    if viz_trees:
        viz_tree(G, 'neural', outdir=figs_dir) 
        viz_tree(sat_tree, 'sat', outdir=figs_dir)
        viz_tree(span_tree, 'mst', outdir=figs_dir)

# ---------------------------------------
    return alphas, mcosts, scosts, first_time

def pareto_analysis(G, neuron_name, neuron_type,\
                    fronts_dir=NEUROMORPHO_FRONTS_DIR,\
                    output_dir=NEUROMORPHO_OUTPUT_DIR,\
                    figs_dir=NEUROMORPHO_FIGS_DIR, output=True,\
                    viz_trees=VIZ_TREES_NEUROMORPHO):
    
    assert G.number_of_nodes() > 0
    assert is_tree(G)

    print neuron_name, neuron_type
   
    #print "making graph"
    synapses = G.graph['synapses']
    points = synapses + [G.graph['root']]
    point_graph = G.subgraph(points)
    print point_graph.number_of_nodes(), 'points'
    point_graph = complete_graph(point_graph)
   
# ---------------------------------------
    tree_costs_fname = '%s/tree_costs.csv' % fronts_dir
   
    models_fname = '%s/models_%s.csv' % (output_dir, neuron_name)
    
    output_fname = '%s/pareto_steiner_%s.csv' %  (output_dir, neuron_name)

# ---------------------------------------
    alphas, mcosts, scosts, first_time = pareto_front(G, point_graph,\
                                                      neuron_name, neuron_type,\
                                                      fronts_dir, figs_dir,\
                                                      viz_trees) 
# ---------------------------------------
    
    neural_mcost = mst_cost(G)
    neural_scost = satellite_cost(G, relevant_nodes=point_graph.nodes())
    
    neural_dist, neural_index = DIST_FUNC(mcosts, scosts, neural_mcost,\
                                          neural_scost)
    neural_closem = mcosts[neural_index] 
    neural_closes = scosts[neural_index]
    neural_alpha = alphas[neural_index]

    
# ---------------------------------------
    centroid_tree = centroid_mst(point_graph)

    centroid_mcost = mst_cost(centroid_tree)
    centroid_scost = satellite_cost(centroid_tree, relevant_nodes=point_graph.nodes())
    
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
        with open(models_fname, 'a') as models_file:
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
        with open(output_fname, 'a') as output_file:
            output_file.write('%s\n' % write_items)

def pareto_analysis_imaris(G, neuron_name, neuron_type,\
                           fronts_dir=IMARIS_FRONTS_DIR,\
                           figs_dir=IMARIS_FIGS_DIR,\
                           viz_trees=VIZ_TREES_IMARIS):
    assert G.number_of_nodes() > 0
    assert is_tree(G)

    print neuron_name
    print "making graph"
    synapses = G.graph['synapses']
    points = synapses + [G.graph['root']]
    point_graph = G.subgraph(points)
    point_graph = complete_graph(point_graph)
    print point_graph.number_of_nodes(), "points"
    
    alphas, mcosts, scosts, first_time = pareto_front(G, point_graph,\
                                                      neuron_name,\
                                                      neuron_type,\
                                                      fronts_dir, figs_dir,\
                                                      viz_trees)
   
    axon = neuron_type == 'axon'
    alphas2, mcosts2, scosts2 = pareto_tree_costs(G, point_graph, axon=axon,\
                                                  viz_trees=False,\
                                                  figs_dir=None)
    mcosts2 = [mcosts[0]] + mcosts2
    scosts2 = [scosts[0]] + scosts2
# ---------------------------------------
    
    neural_mcost = mst_cost(G)
    neural_scost = satellite_cost(G, relevant_nodes=point_graph.nodes())
    
    centroid_tree = centroid_mst(point_graph)
    centroid_mcost = mst_cost(centroid_tree)
    centroid_scost = satellite_cost(centroid_tree, relevant_nodes=point_graph.nodes())
    
    tree_costs_fname = '%s/tree_costs.csv' % fronts_dir
    with open(tree_costs_fname, 'w') as tree_costs_file:
        tree_costs_file.write('tree, mcost, scost\n')
        tree_costs_file.write('%s, %f, %f\n' % ('neural',neural_mcost,\
                                                neural_scost))
        tree_costs_file.write('%s, %f, %f\n' % ('centroid', centroid_mcost,\
                                                centroid_scost))

    pareto_plot(fronts_dir, figs_dir)
    
    pylab.figure()
    sns.set()
    
    pylab.plot(mcosts, scosts, c = 'b', label='_nolegend_')
    pylab.scatter(mcosts, scosts, c='b', label='pareto front')
    
    pylab.plot(mcosts2, scosts2, c = 'm', label='_nolegend_')
    pylab.scatter(mcosts2, scosts2, c='m', label='sandbox pareto front')
    
    pylab.savefig('%s/pareto_front.pdf' % figs_dir, format='pdf')

    neural_dist, neural_index = DIST_FUNC(mcosts, scosts, neural_mcost,\
                                          neural_scost)
    
    centroid_dist, centroid_index = DIST_FUNC(mcosts, scosts, centroid_mcost,\
                                              centroid_scost)
    
    for scale_factor, tree in zip([neural_dist, centroid_dist], ['neural', 'centroid']):
        color = COLORS[tree]
        x = scale_factor * pylab.array(mcosts)
        y = scale_factor * pylab.array(scosts)
        pylab.plot(x, y, c=color, linestyle='-')
        pylab.scatter(x, y, c=color, label='s = %0.2f' % scale_factor)

    pylab.legend()

    pylab.savefig('%s/pareto_front_scaled.pdf' % figs_dir, format='pdf')
    pylab.close()


def pareto_analysis_neuromorpho(min_nodes=MIN_NODES, max_nodes=MAX_NODES,\
                                cell_types=None, animal_species=None,\
                                regions=None, labs=None, names=None,\
                                neuron_types=None, viz_trees=VIZ_TREES_NEUROMORPHO,\
                                plot=False, synthetic=False):
    #directory = 'neuromorpho'
    datasets_dir = NEUROMORPHO_DATASETS_DIR
    for cell_type in os.listdir(datasets_dir):
        if cell_type in SKIP_TYPES:
            continue
        if cell_types != None and cell_type not in cell_types:
            continue
        for species in os.listdir(datasets_dir + '/' + cell_type):
            if animal_species != None and species not in animal_species:
                continue
            for region in os.listdir(datasets_dir + '/' + cell_type + '/' + species):
                if regions != None and region not in regions:
                    continue
                for lab in os.listdir(datasets_dir + "/" + cell_type + '/' + species+ '/' + region):
                    if labs != None and lab not in labs:
                        continue
                    for neuron_file in os.listdir(datasets_dir + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = datasets_dir + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron_file
                        
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

                            if neuron_types != None and i not in neuron_types:
                                continue
                            
                            neuron_type = NEURON_TYPES[i] 
                            
                            H = None
                            if synthetic:
                                H = add_synapses(G, neuron_type=neuron_type)
                            else:
                                H = G.copy()
                                H.graph['synapses'] = []
                                for u in H.nodes_iter():
                                    if u != H.graph['root']:
                                        H.graph['synapses'].append(u)
                         
                            npoints = len(H.graph['synapses']) + 1
                            if not (min_nodes <= npoints <= max_nodes):
                                print "wrong nodes", npoints
                                continue


                            fronts_dir = '%s/%s/%s' % (NEUROMORPHO_FRONTS_DIR, neuron_name, neuron_type)
                            fronts_dir = fronts_dir.replace(' ', '_')

                            output_dir = NEUROMORPHO_TEMP_DIR
                            output_dir = output_dir.replace(' ', '_')

                            figs_dir = '%s/%s/%s/%s/%s/%s/%s' % (NEUROMORPHO_FIGS_DIR,\
                                                                 cell_type,\
                                                                 species,\
                                                                 region, lab,\
                                                                 neuron_name,\
                                                                 neuron_type)

                            figs_dir = figs_dir.replace(' ', '_')
                           
                            if synthetic:
                                fronts_dir = fronts_dir.replace('/pareto_fronts/',\
                                                                '/pareto_fronts_synthetic/')
                                output_dir += '_synthetic'
                                figs_dir = figs_dir.replace('/steiner_figs/',\
                                                            '/steiner_figs_synthetic/')

                            os.system('mkdir -p %s' % fronts_dir)
                            os.system('mkdir -p %s' % output_dir)
                            if plot or viz_trees:
                                os.system('mkdir -p %s' % figs_dir)

                            try: 
                                pareto_analysis(H, neuron_name, neuron_type,\
                                                fronts_dir=fronts_dir,\
                                                output_dir=output_dir,\
                                                figs_dir=figs_dir,\
                                                viz_trees=viz_trees)
                                if plot:
                                    pareto_plot(fronts_dir, figs_dir)


                            except RuntimeError as r:
                                print r
                                continue

def imaris_plots():
    for neuron_name in os.listdir(IMARIS_TRACINGS_DIR):
        if os.path.isdir('%s/%s' % (IMARIS_TRACINGS_DIR, neuron_name)):
            imfiles = []
            for fname in os.listdir('%s/%s' % (IMARIS_TRACINGS_DIR, neuron_name)):
                if 'Position' in fname:
                    imfiles.append('%s/%s/%s' % (IMARIS_TRACINGS_DIR, neuron_name, fname))
            
            G = read_imaris(imfiles, viz=False)
            
            fronts_dir = '%s/%s' % (IMARIS_FRONTS_DIR, neuron_name)
            fronts_dir = fronts_dir.replace(' ', '_')
            os.system('mkdir -p %s' % fronts_dir)

            figs_dir = '%s/%s' % (IMARIS_FIGS_DIR,\
                                  neuron_name)
            figs_dir = figs_dir.replace(' ', '_')
            os.system('mkdir -p %s' % figs_dir)

            neuron_type = None
            if 'axon' in neuron_name:
                neuron_type = 'axon'
            else:
                neuron_type = 'dendrite'


            pareto_analysis_imaris(G, neuron_name, neuron_type,\
                                   fronts_dir=fronts_dir,\
                                   figs_dir=figs_dir,\
                                   viz_trees=VIZ_TREES_IMARIS)

def neuron_builder_analysis(rmin=0.5, rmax=1.5, rstep=0.01, num_iters=10):
    for i in xrange(num_iters):
        for radius in pylab.arange(rmax, rmin - rstep, -rstep):
            print "radius", radius
            G = build_neuron_snider(radius)
            point_graph = complete_graph(G)
            print point_graph.number_of_nodes(), "points"

            print "sorting neighbors"
            sort_neighbors(point_graph)
           
            alphas, mcosts, scosts, = pareto_tree_costs(G, point_graph, axon=False)
            
            neural_mcost = mst_cost(G)
            neural_scost = satellite_cost(G, relevant_nodes=point_graph.nodes())
            
            neural_dist, neural_index = DIST_FUNC(mcosts, scosts, neural_mcost,\
                                                  neural_scost)
            neural_closem = mcosts[neural_index] 
            neural_closes = scosts[neural_index]
            neural_alpha = alphas[neural_index]
            
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
    parser.add_argument('-np', '--neuromorpho_plot', action='store_true')
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
    parser.add_argument('-nt', '--neuron_types', nargs='+', type=int, default=None)
    parser.add_argument('-bt', '--boutons', action='store_true')
    parser.add_argument('--synthetic', action='store_true')

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
    neuron_types = args.neuron_types
    boutons = args.boutons
    synthetic = args.synthetic

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

        pareto_analysis_neuromorpho(G, name, cell_type, species, region, lab)

        return None

    if imaris:
        imaris_plots()
    if neuromorpho or neuromorpho_plot:
        pareto_analysis_neuromorpho(min_nodes, max_nodes, cell_types,\
                                    animal_species, regions, labs, names,\
                                    neuron_types, viz_trees, neuromorpho_plot,\
                                    synthetic)
    if neuron_builder:
        neuron_builder_analysis()
    if boutons:
        boutons_plots()

if __name__ == '__main__':
    main()
