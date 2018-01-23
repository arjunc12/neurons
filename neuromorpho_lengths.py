from neuron_utils import get_neuron_points
import pylab
import os
import argparse
import pandas as pd

def mean_edge_length(G):
    lengths = []
    for u, v in G.edges_iter():
        lengths.append(G[u][v]['length'])
    assert len(lengths) > 0
    return pylab.mean(lengths)

def write_edge_lengths():
    directory = 'datasets'
    i = 0
    f = open('neuromorpho_lengths.csv', 'a')
    for cell_type in os.listdir(directory):
        #if cell_type in SKIP_TYPES:
            #pass
        #if cell_types != None and cell_type not in cell_types:
            #continue
        for species in os.listdir(directory + '/' + cell_type):
            #if animal_species != None and species not in animal_species:
                #continue
            for region in os.listdir(directory + '/' + cell_type + '/' + species):
                #if regions != None and region not in regions:
                    #continue
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
                            if G.number_of_edges() == 0:
                                continue
                            
                            name = neuron[:-8] + str(i)
                            
                            mean_length = mean_edge_length(G)
                            f.write('%s, %d, %f\n' % (name, G.number_of_edges(),\
                                                    mean_length))
    f.close()

def edge_length_stats():
    df = pd.read_csv('neuromorpho_lengths.csv', names=['name', 'points', 'mean_edge_length'])
    print pylab.average(df['mean_edge_length'], weights=df['points'])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--write', action='store_true')
    parser.add_argument('-s', '--stats', action='store_true')

    args = parser.parse_args()
    write = args.write
    stats = args.stats

    if write:
        write_edge_lengths()
    if stats:
        edge_length_stats()

if __name__ == '__main__':
    main()
