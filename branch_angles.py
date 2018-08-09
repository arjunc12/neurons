from neuron_utils import get_neuron_points
import pylab
import os
import pandas as pd
from itertools import combinations

OUTDIR = '/iblsn/data/Arjun/neurons/branch_angles'
FNAME = 'branch_angles.csv'
DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'

NEURON_TYPES = {0 : 'axon', 1 : 'basal dendrite', 2: 'apical dendrite',\
                3: 'truncated axon'}

def branch_angle(G, parent, child1, child2):
    parent_coord = pylab.array(G.node[parent]['coord'])
    child1_coord = pylab.array(G.node[child1]['coord'])
    child2_coord = pylab.array(G.node[child2]['coord'])

    v1 = child1_coord - parent_coord
    v2 = child2_coord - parent_coord

    m1 = ((v1 ** 2).sum()) ** 0.5
    m2 = ((v2 ** 2).sum()) ** 0.5

    dp = pylab.dot(v1, v2)
   
    cos = dp / (m1 * m2)
    if cos < -1:
        cos = 1
    elif cos > 1:
        cos = 1
    theta = pylab.arccos(cos)

    return theta

def get_angles(G):
    root = G.graph['root']
    visited = set()
    queue = [root]
    angles = {}
    while len(queue) > 0:
        curr = queue.pop(0)
        visited.add(curr)
        neighbors = list(G.neighbors(curr))
        children = filter(lambda x : x not in visited, neighbors)
        for child1, child2 in combinations(children, 2):
            angle = branch_angle(G, curr, child1, child2)
            angles[(curr, child1, child2)] = angle
        queue += children

    return angles

def write_angles():
    prev_neurons = set()
    first_line = False
    if FNAME in os.listdir(OUTDIR):
        df = pd.read_csv('%s/%s' % (OUTDIR, FNAME), skipinitialspace=True)
        neuron_names = list(df['neuron name'])
        neuron_types = list(df['neuron type'])
        prev_neurons = set(zip(neuron_names, neuron_types))
    else:
        first_line = True
    i = 0
    with open('%s/%s' % (OUTDIR, FNAME), 'a') as f:
        if first_line:
            f.write('neuron name, neuron type, parent, child1, child2, angle\n')
        directory = DATASETS_DIR
        for cell_type in os.listdir(directory):
            for species in os.listdir(directory + '/' + cell_type):
                for region in os.listdir(directory + '/' + cell_type + '/' + species):
                    for lab in os.listdir(directory + "/" + cell_type + '/' + species+ '/' + region):
                        for neuron in os.listdir(directory + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                            filename = directory + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron
                            
                            if neuron[-8:] != ".CNG.swc": 
                                continue
                           
                            neuron_name = neuron[:-8]

                            try:
                                graphs = get_neuron_points(filename)
                            except AssertionError:
                                continue

                            for i, G in enumerate(graphs):
                                neuron_type = NEURON_TYPES[i]

                                if (neuron_name, neuron_type) in prev_neurons:
                                    continue
                                prev_neurons.add((neuron_name, neuron_type))

                                if G == None:
                                    continue

                                print neuron_name, neuron_type

                                angles = get_angles(G)

                                for (parent, child1, child2), angle in angles.iteritems():
                                    if pylab.isnan(angle):
                                        continue
                                    parent = int(parent)
                                    child1 = int(child1)
                                    child2 = int(child2)
                                    f.write('%s, %s, %d, %d, %d, %f\n' %\
                                            (neuron_name, neuron_type, parent, child1, child2, angle))

def angles_stats():
    df = pd.read_csv('%s/%s' % (OUTDIR, FNAME), skipinitialspace=True)

def main():
    #write_angles()
    angles_stats()

if __name__ == '__main__':
    main()
