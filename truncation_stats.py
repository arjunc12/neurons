from pareto_steiner_stats import get_dfs
import pylab
import os
from neuron_utils import get_neuron_points, viz_tree
import argparse

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'

def truncation_stats():
    df = get_df()
    df = df.drop_duplicates(subset='name')
    df['swc'] = df['name'].str.slice(stop=-1)
    #df = df[['swc', 'neuron_type', 'points']]
    original_points = []
    new_points = []
    for name, group in df.groupby('swc'):
        points = group['points']
        neuron_types = group['neuron_type']
        index1 = neuron_types == 'axon'
        index2 = neuron_types == 'truncated axon'
        if 'axon' in neuron_types.values and 'truncated axon' in neuron_types.values:
            points1 = list(points[index1])[0]
            points2 = list(points[index2])[0]
            original_points.append(float(points1))
            new_points.append(float(points2))

    original_points = pylab.array(original_points)
    new_points = pylab.array(new_points)
    diff = original_points - new_points
    print pylab.mean(diff)
    print pylab.mean(diff / original_points)

def truncation_plots():
    directory = DATASETS_DIR
    i = 0
    for cell_type in os.listdir(directory):
        for species in os.listdir(directory + '/' + cell_type):
            for region in os.listdir(directory + '/' + cell_type + '/' + species):
                for lab in os.listdir(directory + "/" + cell_type + '/' + species+ '/' + region):
                    for neuron in os.listdir(directory + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = directory + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron
                        
                        if neuron[-8:] != ".CNG.swc": 
                            continue
                        name = neuron[:-8]
                        
                        graphs = None
                        try:
                            graphs = get_neuron_points(filename)
                        except AssertionError as e:
                            continue

                        if len(graphs) != 4:
                            print graphs
                        assert len(graphs) == 4
                        axon, truncated_axon = graphs[0], graphs[3]
                        if axon != None and truncated_axon != None:
                            if axon.number_of_nodes() != truncated_axon.number_of_nodes():
                                if truncated_axon.number_of_nodes() > 1:
                                    if i <= 20:
                                        print i
                                        outdir = 'truncation/drawings'
                                        name1 = name + '-axon1'
                                        name2 = name + '-axon2'
                                        print name1, name2, outdir
                                        '''
                                        axon.graph['synapses'] = []
                                        truncated_axon.graph['synapses'] = []
                                        for u in axon.nodes():
                                            axon.graph['synapses'].append(u)
                                            axon.node[u]['label'] = 'synapse'
                                            if truncated_axon.has_node(u):
                                                truncated_axon.graph['synapses'].append(u)
                                                truncated_axon.node[u]['label'] = 'synapse'
                                        '''
                                        root1 = axon.graph['root']
                                        root2 = truncated_axon.graph['root']
                                        axon.node[root1]['label'] = 'root'
                                        truncated_axon.node[root2]['label'] = 'root'
                                        kwargs = viz_tree(axon, name1, outdir)
                                        viz_tree(truncated_axon, name2, outdir, **kwargs)
                                        i += 1
                                        #exit()
                                    else:
                                        exit()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--stats', action='store_true')
    parser.add_argument('-p', '--plots', action='store_true')

    args = parser.parse_args()
    stats = args.stats
    plots = args.plots

    if stats:
        truncation_stats()
    if plots:
        truncation_plots()

if __name__ == '__main__':
    main()
