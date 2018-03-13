from neuron_utils import get_neuron_points, add_synapses
import pandas as pd
import argparse
import os

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'
OUTDIR = '/iblsn/data/Arjun/neurons/neuron_size'

NEURON_TYPES = {0 : 'axon', 1 : 'basal dendrite', 2: 'apical dendrite',\
                3: 'truncated axon'}

def get_sizes(synthetic=False):
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--synthetic', action='store_true')
    args = parser.parse_args()
    synthetic = args.synthetic

    get_sizes(synthetic=synthetic)

if __name__ == '__main__':
    main()
