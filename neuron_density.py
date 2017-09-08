import matplotlib as mpl
mpl.use('agg')
import pylab
from cost_functions import mst_cost
from scipy.spatial import ConvexHull
from scipy.spatial.qhull import QhullError
from neuron_utils import get_neuron_points
import pandas as pd
import argparse
import os
from pareto_mst import SKIP_TYPES

def graph_volume(G):
    points = []
    for u in G.nodes_iter():
        points.append(G.node[u]['coord'])
    hull = ConvexHull(points)
    return hull.volume

def density(G):
    return mst_cost(G) / graph_volume(G)

def get_densities():
    prev_names = set()
    if 'neuron_density.csv' in os.listdir('.'):
        df = pd.read_csv('neuron_density.csv', names=['name', 'points', 'volume', 'mcost'])
        prev_names = set(df['name'].values)
    directory = 'datasets'
    i = 0
    f = open('neuron_density.csv', 'a')
    for cell_type in os.listdir(directory):
        if cell_type in SKIP_TYPES:
            continue
        for species in os.listdir(directory + '/' + cell_type):
            for region in os.listdir(directory + '/' + cell_type + '/' + species):
                for lab in os.listdir(directory + "/" + cell_type + '/' + species+ '/' + region):
                    for neuron in os.listdir(directory + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = directory + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron
                        
                        if neuron[-8:] != ".CNG.swc": 
                            continue
                        
                        try:
                            graphs = get_neuron_points(filename)
                        except AssertionError:
                            continue

                        for i, G in enumerate(graphs):
                            name = neuron[:-8] + str(i)

                            if name in prev_names:
                                continue

                            if G.number_of_nodes() < 4:
                                continue

                            print name
                            volume = None
                            try:
                                volume = graph_volume(G)
                            except QhullError as qerror:
                                continue
                            
                            if volume == None:
                                continue

                            mcost = mst_cost(G)

                            write_items = [name, G.number_of_nodes(), volume, mcost]
                            write_items = map(str, write_items)
                            write_items = ', '.join(write_items)

                            f.write('%s\n' % write_items)

    f.close()

def get_df():
    df = pd.read_csv('neuron_density.csv',\
                      names=['name', 'points', 'volume', 'mcost'],\
                      skipinitialspace=True)
    return df

def neuron_density_stats():
    df = get_df()
    pylab.figure()
    #pylab.scatter(df['volume'], df['mcost'])
    pylab.scatter(pylab.log(df['volume']), pylab.log(df['mcost']))
    pylab.xlabel('volume')
    pylab.ylabel('strategy')
    pylab.savefig('neuron_density.pdf', format='pdf')
    pylab.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--density', action='store_true')
    parser.add_argument('-s', '--stats', action='store_true')
    
    args = parser.parse_args()
    density = args.density
    stats = args.stats

    if stats:
        neuron_density_stats()

    if density:
        get_densities()

if __name__ == '__main__':
    main()
