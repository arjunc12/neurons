import os
import matplotlib as mpl
mpl.use('agg')
import networkx as nx
from neuron_utils import *
from collections import defaultdict
import numpy as np
import pylab
from graph_utils import is_tree
import argparse

MAX_SEGMENTS = float("inf")

VIZ_TREE = False

def read_imaris(trace_pos, viz=VIZ_TREE, outname='imaris'):
    nodeid = 0
    G = nx.Graph()
    G.graph['synapses'] = []
    start_counts = defaultdict(int)
    coord_ids = {}
    
    start_ids = []
    end_ids = []

    start_end = {}

    for i, pos in enumerate(trace_pos):
        if i + 1 > MAX_SEGMENTS:
            break
        
        num_lines = num_lines = sum(1 for line in open(pos)) - 1
        pos_file = open(pos)
        coords = []
        k = 0
        for row in pos_file:
            if k == 0:
                k += 1
                continue
            row = row.split(',')
            x = float(row[0])
            y = float(row[1])
            z = float(row[2])
            coord = (x, y, z)
            coord_id = int(row[-2])
            coords.append((coord_id, coord))

        coords = sorted(coords)
        prev_coord = None
        for j, (coord_id, new_coord) in enumerate(coords):
            id1 = None
            if new_coord in coord_ids:
                id1 = coord_ids[new_coord]
            else:
                id1 = nodeid
                coord_ids[new_coord] = id1
            G.add_node(id1)
            G.node[id1]['coord'] = new_coord
            G.node[id1]['label'] = 'synapse'
            G.graph['synapses'].append(id1)
            if j == 0:
                start_counts[new_coord] += 1
                start_ids.append(id1)
            else:
                assert prev_coord != None
                id2 = coord_ids[prev_coord]
                if id1 != id2:
                    G.add_edge(id1, id2)
                    G[id1][id2]['length'] = point_dist(G.node[id1]['coord'], G.node[id2]['coord'])

                if j == len(coords) - 1:
                    end_ids.append(id1)
            
            nodeid += 1
            prev_coord = new_coord

        pos_file.close()

    root_coord =  max(start_counts.keys(), key = (lambda x : start_counts[x]))
    root_id = coord_ids[root_coord]
    G.graph['root'] = root_id
    G.node[root_id]['label'] = 'root'

    for start_id in start_ids:
        if start_id != root_id:
            closest_dist = float("inf")
            closest_end = None
            for end_id in end_ids:
                if nx.has_path(G, start_id, end_id):
                    continue
                start_coord = G.node[start_id]['coord']
                end_coord = G.node[end_id]['coord']
                dist = point_dist(start_coord, end_coord)
                if dist < closest_dist:
                    closest_dist = dist
                    closest_end = end_id
            G.add_edge(start_id, closest_end)
            G[start_id][closest_end]['length'] = closest_dist

    connected_components = nx.connected_components(G)
    connected_components = list(connected_components)
    isolated_nodes = []
    for component in connected_components:
        if G.graph['root'] not in component:
            isolated_nodes += component


    for u, v in G.edges():
        if u == v:
            G.remove_edge(u, v)

    label_points(G)

    for start_id in start_ids:
        if start_id in isolated_nodes:
            G.node[start_id]['label'] = 'isolated_start'

    for end_id in end_ids:
        if end_id in isolated_nodes:
            G.node[end_id]['label'] = 'isolated_end'

    if viz:
        viz_tree(G, outname, 'imaris')
    print G.number_of_nodes()
    print G.number_of_edges()
    print nx.is_connected(G)

    f = open('imaris_lengths.csv', 'a')
    edge_lengths = []
    for u, v in G.edges_iter():
        length = G[u][v]['length']
        edge_lengths.append(length)
        f.write('%s, %f\n' % (outname, length))
    f.close()
    print "mean edge length", pylab.mean(edge_lengths)

    return G

def draw_imaris():
    for i, neuron in enumerate(os.listdir('imaris')):
        print neuron
        if os.path.isdir('imaris/%s' % neuron):
            trace_pos = []
            for fname in os.listdir('imaris/%s' % neuron):
                if 'Position' in fname:
                    trace_pos.append('imaris/%s/%s' % (neuron, fname))
            read_imaris(trace_pos, viz=True, outname=neuron)

def plot_edge_lengths():
    import pandas as pd
    df = pd.read_csv('imaris_lengths.csv', names = ['neuron', 'length'])
    pylab.hist(df['length'])
    pylab.savefig('imaris/imaris_lengths.pdf', format='pdf')
    pylab.close()
    mu = pylab.mean(df['length'])
    sigma2 = pylab.var(df['length'], ddof=1)
    l = 1.0 / mu
    l2 = l ** 2
    print "mean", mu
    print "lambda", l
    print 'variance', sigma2
    print 'var-hat', 1.0 / l2

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--draw', action='store_true')
    parser.add_argument('-l', '--lengths', action='store_true')

    args = parser.parse_args()
    draw = args.draw
    lengths = args.lengths

    if draw or not lengths:
        draw_imaris()
    else:
        plot_edge_lengths()

if __name__ == '__main__':
    main()
