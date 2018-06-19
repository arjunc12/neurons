import pandas as pd
import os
from dist_functions import pareto_dist_scale

def get_neural_costs(trees_df):
    trees_df = trees_df[trees_df['tree'] == 'neural']
    mcost = list(trees_df['mcost'])[0]
    scost = list(trees_df['scost'])[0]
    return mcost, scost

fronts_dir = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts_synthetic'
tree_file = 'tree_costs.csv'

with open('/iblsn/data/Arjun/neurons/pareto_steiner_output/neural_costs.csv', 'w') as models_file:
    models_file.write('neuron name, neuron type, model, dist, success, ratio\n')
    for neuron_name in os.listdir(fronts_dir):
        for neuron_type in os.listdir('%s/%s' % (fronts_dir, neuron_name)):
            files = os.listdir('%s/%s/%s' % (fronts_dir, neuron_name, neuron_type))

            if tree_file not in files:
                continue
            trees_df = pd.read_csv('%s/%s' % (front_dir, tree_file), skipinitialspace=True)
            
            print neuron_name, neuron_type
            
            neural_mcost, neural_scost = get_neural_costs(trees_df)

