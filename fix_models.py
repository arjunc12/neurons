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
parser.add_argument('--synthetic', action='store_true')

args = parser.parse_args()
synthetic = args.synthetic

pareto_fronts_dir = NEUROMORPHO_FRONTS_DIR
output_dir = NEUROMORPHO_TEMP_DIR
if synthetic:
    pareto_fronts_dir += '_synthetic'
    output_dir += '_synthetic'

print pareto_fronts_dir
for neuron_name in os.listdir(pareto_fronts_dir):
    for neuron_type in os.listdir('%s/%s' % (pareto_fronts_dir, neuron_name)):
        fronts_dir = '%s/%s/%s' % (pareto_fronts_dir, neuron_name, neuron_type)

        print neuron_name, neuron_type
        ntype = neuron_type.replace('_', ' ')
     
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
        
        trees_df = pd.read_csv('%s/tree_costs.csv' % fronts_dir, skipinitialspace=True)
        neural_tree = trees_df[trees_df['tree'] == 'neural']
        neural_mcost = list(neural_tree['mcost'])[0]
        neural_scost = list(neural_tree['scost'])[0]
        print "neural costs", neural_mcost, neural_scost
        
        neural_dist, neural_index = DIST_FUNC(mcosts, scosts, neural_mcost,\
                                              neural_scost)
        neural_closem = mcosts[neural_index] 
        neural_closes = scosts[neural_index]
        neural_alpha = alphas[neural_index]

        
        with open(models_fname, 'a') as models_file:
            models_file.write('%s, %s, %s, %f\n' % (neuron_name,\
                                                      neuron_type,\
                                                      'neural',\
                                                      neural_dist))
            
            for tree, group in trees_df.groupby('tree'):
                for null_mcost, null_scost in zip(group['mcost'], group['scost']):
                    null_dist, null_index = DIST_FUNC(mcosts, scosts, null_mcost, null_scost)
                    null_success = null_dist <= neural_dist
                    null_ratio = null_dist / neural_dist
                    models_file.write('%s, %s, %s, %f, %d, %f\n' % (neuron_name,\
                                      neuron_type, tree, null_dist,\
                                      null_success, null_ratio))
