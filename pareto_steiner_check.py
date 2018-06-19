import pandas as pd
import os
from dist_functions import pareto_dist_scale
from make_neuron_pbs import read_sizes

def get_neural_costs(trees_df):
    trees_df = trees_df[trees_df['tree'] == 'neural']
    mcost = list(trees_df['mcost'])[0]
    scost = list(trees_df['scost'])[0]
    return mcost, scost

sizes = read_sizes(size_file='/iblsn/data/Arjun/neurons/neuron_size/neuron_size_synthetic.csv')
fronts_dir = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts_synthetic'
tree_file = 'tree_costs.csv'
front_file = 'pareto_front.csv'

no_trees = open('no_trees.csv', 'w')
no_size = open('no_size.csv', 'w')
with open('/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_steiner_synthetic_fixed.csv', 'w') as output_file:
    output_file.write('neuron name, neuron type, points, alpha\n')
    for neuron_name in os.listdir(fronts_dir):
        for neuron_type in os.listdir('%s/%s' % (fronts_dir, neuron_name)):
            ntype = neuron_type.replace('_', ' ')
            if (neuron_name, ntype) not in sizes:
                no_size.write('%s, %s\n' % (neuron_name, ntype))
                print neuron_name, ntype, "no size"
                continue
            size = sizes[(neuron_name, ntype)]
            files = os.listdir('%s/%s/%s' % (fronts_dir, neuron_name, neuron_type))
            if front_file not in files:
                continue

            front_dir = '%s/%s/%s' % (fronts_dir, neuron_name, neuron_type)
            front_df = pd.read_csv('%s/%s' % (front_dir, front_file), skipinitialspace=True)
            opt_mcost = front_df['mcost']
            opt_scost = front_df['scost']
            alpha = front_df['alpha']

            if tree_file not in files:
                no_trees.write('%s, %s\n' % (neuron_name, ntype))
                continue
            trees_df = pd.read_csv('%s/%s' % (front_dir, tree_file), skipinitialspace=True)
            
            print neuron_name, ntype
            
            neural_mcost, neural_scost = get_neural_costs(trees_df)
            neural_dist, neural_index = pareto_dist_scale(opt_mcost, opt_scost, neural_mcost, neural_scost)
            neural_alpha = alpha[neural_index]
            output_file.write('%s, %s, %d, %f\n' % (neuron_name, ntype, size, neural_alpha))

no_trees.close()
