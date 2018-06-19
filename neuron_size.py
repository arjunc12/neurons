from neuron_utils import get_neuron_points, add_synapses
import pandas as pd
import argparse
import os

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'
OUTDIR = '/iblsn/data/Arjun/neurons/neuron_size'
FRONTS_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts'
SIZE_FILE = '/iblsn/data/Arjun/neurons/neuron_size/neuron_size.csv'
MIN_POINTS = 50

NEURON_TYPES = {0 : 'axon', 1 : 'basal dendrite', 2: 'apical dendrite',\
                3: 'truncated axon'}

def write_points(synthetic=False):
    prev_neurons = set()
    first_line = False
    fname = 'neuron_size'
    if synthetic:
        fname += '_synthetic'
    fname += '.csv'
    if fname in os.listdir(OUTDIR):
        df = pd.read_csv('%s/%s' % (OUTDIR, fname),\
                          names=['neuron name', 'neuron type', 'points'],\
                          skipinitialspace=True)
        neuron_names = list(df['neuron name'])
        neuron_types = list(df['neuron type'])
        prev_neurons = set(zip(neuron_names, neuron_types))
    else:
        first_line = True
    i = 0
    with open('%s/%s' % (OUTDIR, fname), 'a') as f:
        if first_line:
            f.write('neuron name, neuron type, points\n')
        directory = DATASETS_DIR
        for cell_type in os.listdir(directory):
            for species in os.listdir(directory + '/' + cell_type):
                for region in os.listdir(directory + '/' + cell_type + '/' + species):
                    for lab in os.listdir(directory + "/" + cell_type + '/' + species+ '/' + region):
                        for neuron in os.listdir(directory + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                            filename = directory + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron
                            if filename[-8:] != '.CNG.swc':
                                continue
                           
                            neuron_name = neuron[:-8]

                            try:
                                graphs = get_neuron_points(filename)
                            except AssertionError:
                                print "exception"
                                continue
                            print "no exception"

                            for i, G in enumerate(graphs):
                                neuron_type = NEURON_TYPES[i]
                                
                                print neuron_name, neuron_type
                                print len(prev_neurons)

                                if (neuron_name, neuron_type) in prev_neurons:
                                    continue

                                if G == None:
                                    continue
                                
                                prev_neurons.add((neuron_name, neuron_type))

                                if synthetic:
                                    G = add_synapses(G, neuron_type=neuron_type)
                                else:
                                    G.graph['synapses'] = []
                                    for u in G.nodes():
                                        if u != G.graph['root']:
                                            G.graph['synapses'].append(u)

                                points = len(G.graph['synapses']) + 1
                                write_items = [neuron_name, neuron_type, points]
                                
                                write_items = map(str, write_items)
                                write_items = ', '.join(write_items)

                                f.write('%s\n' % write_items)

def write_neural_costs(synthetic=False):
    fronts_dir = FRONTS_DIR
    if synthetic:
        fronts_dir += '_synthetic'
    
    costs_fname = '/iblsn/data/Arjun/neurons/pareto_steiner_output/neural_costs.csv'
    if synthetic:
        costs_fname = costs_fname.replace('.csv', '_synthetic.csv')
    print costs_fname
    with open(costs_fname, 'w') as costs_file:
        costs_file.write('neuron name, neuron type, mcost, scost\n')
        for neuron_name in os.listdir(fronts_dir):
            for neuron_type in os.listdir('%s/%s' % (fronts_dir, neuron_name)):
                front_dir = '%s/%s/%s' % (fronts_dir, neuron_name, neuron_type)
                if 'tree_costs.csv' not in os.listdir(front_dir):
                    continue
                print neuron_name, neuron_type
                tree_costs = pd.read_csv('%s/tree_costs.csv' % front_dir, skipinitialspace=True)
                tree_costs = tree_costs[tree_costs['tree'] == 'neural']
                mcost = list(tree_costs['mcost'])[0]
                scost = list(tree_costs['scost'])[0]
                costs_file.write('%s, %s, %f, %f\n' % (neuron_name, neuron_type, mcost, scost))

    '''
    costs_fname = '/iblsn/data/Arjun/neurons/pareto_steiner_output/neural_costs.csv'
    if synthetic:
        costs_fname = costs_fname.replace('.csv', '_synthetic.csv')
    costs_df.to_csv(costs_fname, index=False)

    size_file = SIZE_FILE
    if synthetic:
        size_file = size_file.replace('.csv', '_synthetic.csv')
    size_df = pd.read_csv(size_file)
    size_df = size_df[size_df['points'] >= MIN_POINTS]
    
    costs_df = costs_df.merge(size_df)
    print "mcost range", min(costs_df['mcost']), max(costs_df['mcost'])
    print "scost range", min(costs_df['scost']), max(costs_df['scost'])
    '''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--synthetic', action='store_true')
    parser.add_argument('--points', action='store_true')
    parser.add_argument('--costs', action='store_true')
    
    args = parser.parse_args()
    synthetic = args.synthetic
    points = args.points
    costs = args.costs

    if points:
        write_points(synthetic=synthetic)
    if costs:
        write_neural_costs(synthetic=synthetic)

if __name__ == '__main__':
    main()
