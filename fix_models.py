import argparse
import os
from neuron_utils import get_neuron_points, add_synapses
from cost_functions import mst_cost, satellite_cost, graph_costs
from graph_utils import is_tree
from dist_functions import pareto_dist_scale
import pandas as pd
import pylab
from pareto_functions import centroid_mst

SKIP_TYPES = ['Unknown_neurotransmitter', 'Not_reported']

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

NEURON_BUILDER_TREE_DIR = '%s/neuron_builder' % DATA_DRIVE
SNIDER_DIR = '%s/snider' % NEURON_BUILDER_TREE_DIR
SNIDER_FIGS_DIR = '%s/figs' % SNIDER_DIR
SNIDER_TREES_DIR = '%s/trees'

NEURON_TYPES = {0 : 'axon', 1 : 'basal dendrite', 2: 'apical dendrite',\
                3: 'truncated axon'}

ARBOR_TYPES = {'axon': 0, 'basal_dendrite' : 1, 'apical_dendrite' : 2, 'truncated_axon' : 3}

DIST_FUNC = pareto_dist_scale

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cell_types', nargs='+', default=None)
parser.add_argument('-s', '-a', '--species', '--animal_species,', nargs='+',\
                    default=None, dest='animal_species')
parser.add_argument('-l', '--labs', nargs='+', default=None, dest='labs')
parser.add_argument('-na', '--names', nargs='+', default=None, dest='names')
parser.add_argument('-r', '--regions', nargs='+', default=None)
parser.add_argument('-nt', '--neuron_types', nargs='+', default=None)
parser.add_argument('--synthetic', action='store_true')

args = parser.parse_args()
cell_types = args.cell_types
animal_species = args.animal_species
regions = args.regions
labs = args.labs
names = args.names
ntypes = args.neuron_types
neuron_types = None
if ntypes != None:
    neuron_types = []
    for ntype in ntypes:
        if ntype.isdigit():
            neuron_types.append(int(ntype))
        elif ntype in ARBOR_TYPES:
            neuron_types.append(ARBOR_TYPES[ntype])
        else:
            raise ValueError('invalid neuron type')
synthetic = args.synthetic

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
                        
                        fronts_dir = '%s/%s/%s' % (NEUROMORPHO_FRONTS_DIR, neuron_name, neuron_type)
                        fronts_dir = fronts_dir.replace(' ', '_')

                        output_dir = NEUROMORPHO_TEMP_DIR
                        output_dir = output_dir.replace(' ', '_')

                        if synthetic:
                            fronts_dir = fronts_dir.replace('/pareto_fronts/',\
                                                            '/pareto_fronts_synthetic/')
                            output_dir += '_synthetic'
                                
                        assert H.number_of_nodes() > 0
                        assert is_tree(H)

                        print neuron_name, neuron_type
                       
                        synapses = H.graph['synapses']
                        points = synapses + [H.graph['root']]
                        point_graph = H.subgraph(points)
                        print point_graph.number_of_nodes(), 'points'
                        
                        tree_costs_fname = '%s/tree_costs.csv' % fronts_dir

                        models_fname = '%s/models_%s.csv' % (output_dir, neuron_name)
                        
                        pareto_front_fname = '%s/pareto_front.csv' % fronts_dir

                        alphas = None
                        mcosts = None
                        scosts = None

                        df = pd.read_csv('%s/pareto_front.csv' % fronts_dir, skipinitialspace=True)
                        alphas = pylab.array(df['alpha'])
                        mcosts = pylab.array(df['mcost'])
                        scosts = pylab.array(df['scost'])
                         
                        #neural_mcost = mst_cost(H)
                        #neural_scost = satellite_cost(H, relevant_nodes=point_graph.nodes())
                        neural_mcost, neural_scost = graph_costs(H, relevant_nodes=point_graph.nodes())
                        print "neural costs", neural_mcost, neural_scost
                        
                        neural_dist, neural_index = DIST_FUNC(mcosts, scosts, neural_mcost,\
                                                              neural_scost)
                        neural_closem = mcosts[neural_index] 
                        neural_closes = scosts[neural_index]
                        neural_alpha = alphas[neural_index]

                    # ---------------------------------------
                        centroid_tree = centroid_mst(point_graph)

                        #centroid_mcost = mst_cost(centroid_tree)
                        #centroid_scost = satellite_cost(centroid_tree, relevant_nodes=point_graph.nodes())
                        centroid_mcost, centroid_scost = graph_costs(centroid_tree, relevant_nodes=point_graph.nodes())
                        print "centroid costs", centroid_mcost, centroid_scost
                        
                        centroid_dist, centroid_index = DIST_FUNC(mcosts, scosts,\
                                                                  centroid_mcost,\
                                                                  centroid_scost)
                        centroid_closem = mcosts[centroid_index]
                        centroid_closes = scosts[centroid_index]
                        centroid_alpha = alphas[centroid_index]


                    # ---------------------------------------
                        centroid_success = int(centroid_dist <= neural_dist)
                        centroid_ratio = centroid_dist / neural_dist

                        print neuron_name, neuron_type
                        print "neural", neural_dist
                        print "centroid", centroid_dist, centroid_success, centroid_ratio
                        exit()
                        
                        with open(models_fname, 'a') as models_file:
                            models_file.write('%s, %s, %s, %f\n' % (neuron_name,\
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
                        with open(tree_costs_fname, 'w') as tree_costs_file:
                            tree_costs_file.write('tree, mcost, scost\n')
                            tree_costs_file.write('%s, %f, %f\n' % ('neural', neural_mcost,\
                                                                    neural_scost))
                            tree_costs_file.write('%s, %f, %f\n' % ('centroid', centroid_mcost,\
                                                                    centroid_scost))
