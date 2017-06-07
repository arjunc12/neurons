import os
import matplotlib as mpl
mpl.use('agg')
import pandas as pd
import networkx as nx
from neuron_utils import *

def read_imaris(trace_pos):
    nodeid = 0
    start_points = []
    G = nx.Graph()
    for i, pos in enumerate(trace_pos):
        df = pd.read_csv(pos)
        df.sort_values(by='ID', inplace=True)
        rows = list(df.iterrows())
        
        for j in xrange(len(rows) - 1):
            index1, row1 = rows[j]
            index2, row2 = rows[j + 1]
            id1, id2 = nodeid, nodeid + 1
            if j == 0:
                if i == 0:
                    G.add_node(id1)
                    G.node[id1]['coord'] = (row1['Position X'], row1['Position Y'], row1['Position Z'])
                    G.graph['root'] = id1
                else:
                    id1 = G.graph['root']
            G.add_edge(id1, id2)
            G.node[id2]['coord'] = (row2['Position X'], row2['Position Y'], row2['Position Z'])
            G[id1][id2]['length'] = point_dist(G.node[id1]['coord'], G.node[id2]['coord'])
            nodeid += 1

    label_points(G)
    #viz_tree(G, 'imaris', 'imaris')
    return G

def main():
    trace_pos = []
    for fname in os.listdir('imaris'):
        if 'Position' in fname:
            trace_pos.append('imaris/' + fname)
    read_imaris(trace_pos)

if __name__ == '__main__':
    main()
