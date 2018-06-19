import pandas as pd
import os
from dist_functions import pareto_dist_scale

'''
neuron_name = '027-green1'
neuron_type = 'axon'

front_dir = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts_synthetic/%s/%s' % (neuron_name, neuron_type)
trees_df = pd.read_csv('%s/tree_costs.csv' % front_dir, skipinitialspace=True)
front_df = pd.read_csv('%s/pareto_front.csv' % front_dir, skipinitialspace=True)

opt_mcost = front_df['mcost']
opt_scost = front_df['scost']

for tree, group in trees_df.groupby('tree'):
    print tree
    mcosts, scosts = group['mcost'], group['scost']
    for mcost, scost in zip(mcosts, scosts):
        print pareto_dist_scale(opt_mcost, opt_scost, mcost, scost)
'''
def get_neural_costs(trees_df):
    trees_df = trees_df[trees_df['tree'] == 'neural']
    mcost = list(trees_df['mcost'])[0]
    scost = list(trees_df['scost'])[0]
    return mcost, scost

fronts_dir = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts_synthetic'
tree_file = 'tree_costs.csv'
front_file = 'pareto_front.csv'

no_trees = open('no_trees.csv', 'w')
few_trees = open('few_trees.csv', 'w')
with open('/iblsn/data/Arjun/neurons/pareto_steiner_output/models_synthetic_fixed.csv', 'w') as models_file:
    models_file.write('neuron name, neuron type, model, dist, success, ratio\n')
    for neuron_name in os.listdir(fronts_dir):
    #for neuron_name in ['027-green1']:
        for neuron_type in os.listdir('%s/%s' % (fronts_dir, neuron_name)):
        #for neuron_type in ['axon']:
            files = os.listdir('%s/%s/%s' % (fronts_dir, neuron_name, neuron_type))
            if front_file not in files:
                continue

            front_dir = '%s/%s/%s' % (fronts_dir, neuron_name, neuron_type)
            front_df = pd.read_csv('%s/%s' % (front_dir, front_file), skipinitialspace=True)
            opt_mcost = front_df['mcost']
            opt_scost = front_df['scost']

            if tree_file not in files:
                no_trees.write('%s, %s\n' % (neuron_name, neuron_type))
                continue
            trees_df = pd.read_csv('%s/%s' % (front_dir, tree_file), skipinitialspace=True)
            if len(trees_df['tree'].unique()) < 4:
                few_trees.write('%s, %s\n' % (neuron_name, neuron_type))
                continue
            
            print neuron_name, neuron_type
            
            neural_mcost, neural_scost = get_neural_costs(trees_df)
            neural_dist, neural_index = pareto_dist_scale(opt_mcost, opt_scost, neural_mcost, neural_scost)
            ntype = neuron_type.replace('_', ' ')
            models_file.write('%s, %s, %s, %f,,\n' % (neuron_name, ntype, 'neural', neural_dist))

            for tree, group in trees_df.groupby('tree'):
                if tree == 'neural':
                    continue
                for mcost, scost in zip(group['mcost'], group['scost']):
                    dist, index = pareto_dist_scale(opt_mcost, opt_scost, mcost, scost)
                    success = int(dist <= neural_dist)
                    ratio = dist / neural_dist
                    models_file.write('%s, %s, %s, %f, %d, %f\n' % (neuron_name, ntype, tree, dist, success, ratio))
no_trees.close()
few_trees.close()
