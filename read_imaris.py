import os
import matplotlib as mpl
mpl.use('agg')
import pandas as pd
import networkx as nx
from neuron_utils import *
from collections import defaultdict

MAX_SEGMENTS = 30

def read_imaris(trace_pos):
    nodeid = 0
    G = nx.Graph()
    start_counts = defaultdict(int)
    coord_ids = {}
    
    start_ids = []
    end_ids = []

    for i, pos in enumerate(trace_pos):
        if i + 1 > MAX_SEGMENTS:
            break
        df = pd.read_csv(pos)
        df.sort_values(by='ID', inplace=True)
        rows = list(df.iterrows())
        for j in xrange(len(rows)):
            index1, row1 = rows[j]
            id1 = None
            new_coord = (row1['Position X'], row1['Position Y'], row1['Position Z'])
            if new_coord in coord_ids:
                id1 = coord_ids[new_coord]
            else:
                id1 = nodeid
                coord_ids[new_coord] = id1
            G.add_node(id1)
            G.node[id1]['coord'] = new_coord
            if j == 0:
                start_counts[new_coord] += 1
                start_ids.append(id1)
            else:
                index2, row2 = rows[j - 1]
                prev_coord = (row2['Position X'], row2['Position Y'], row2['Position Z'])
                id2 = coord_ids[prev_coord]
                G.add_edge(id1, id2)
                G[id1][id2]['length'] = point_dist(G.node[id1]['coord'], G.node[id2]['coord'])

                if j == len(rows) - 1:
                    end_ids.append(id1)
            
            nodeid += 1

    root_coord =  max(start_counts.keys(), key = (lambda x : start_counts[x]))
    root_id = coord_ids[root_coord]
    G.graph['root'] = root_id

    for start_id in start_ids:
        if start_id != root_id:
            closest_dist = float("inf")
            closest_end = None
            for end_id in end_ids:
                start_coord = G.node[start_id]['coord']
                end_coord = G.node[end_id]['coord']
                dist = point_dist(start_coord, end_coord)
                if dist < closest_dist:
                    closest_dist = dist
                    closest_end = end_id
            G.add_edge(start_id, closest_end)
            G[start_id][closest_end]['length'] = closest_dist

    label_points(G)
    #viz_tree(G, 'imaris2', 'imaris2')
    #print G.number_of_nodes()
    #print G.number_of_edges()
    return G

def main():
    trace_pos = []
    for fname in os.listdir('imaris2'):
        if 'Position' in fname:
            trace_pos.append('imaris2/' + fname)
    read_imaris(trace_pos)

if __name__ == '__main__':
    main()
