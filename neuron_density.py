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

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'
OUTDIR = '/iblsn/data/Arjun/neurons/neuron_density'

NEURON_TYPES = {0 : 'axon', 1 : 'basal dendrite', 2: 'apical dendrite',\
                3: 'truncated axon'}

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
    first_line = False
    if 'neuron_density.csv' in os.listdir(OUTDIR):
        df = pd.read_csv('%s/neuron_density.csv' % OUTDIR,\
                          names=['name', 'points', 'volume', 'mcost'],\
                          skipinitialspace=True)
        prev_names = set(df['name'].values)
    else:
        first_line = True
    i = 0
    with open('%s/neuron_density.csv' % OUTDIR, 'a') as f:
        if first_line:
            f.write('neuron name, neuron type, points, mcost, volume\n')
        directory = DATASETS_DIR
        for cell_type in os.listdir(directory):
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
                                neuron_name = neuron[:-8] + str(i)
                                neuron_type = NEURON_TYPES[i]

                                if neuron_name in prev_names:
                                    continue

                                if G == None:
                                    continue

                                if G.number_of_nodes() < 4:
                                    continue

                                print neuron_name
                                volume = None
                                try:
                                    volume = graph_volume(G)
                                except QhullError as qerror:
                                    continue
                                
                                if volume == None:
                                    continue

                                mcost = mst_cost(G)

                                write_items = [neuron_name, neuron_type, G.number_of_nodes(), mcost, volume]
                                write_items = map(str, write_items)
                                write_items = ', '.join(write_items)

                                f.write('%s\n' % write_items)

def get_df():
    df = pd.read_csv('%s/neuron_density.csv' % OUTDIR,\
                      skipinitialspace=True)
    df['density'] = df['mcost'] / df['volume']
    return df

def neuron_density_stats():
    df = get_df()
    x = pylab.log10(pylab.array(df['volume']))
    y = pylab.log10(pylab.array(df['density']))
    pylab.figure()
    pylab.scatter(x, y)
    pylab.xlabel('log10-volume')
    pylab.ylabel('log10-density')
    pylab.savefig('neuron_density/neuron_density.pdf', format='pdf')
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
