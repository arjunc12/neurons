import matplotlib as mpl
mpl.use('agg')
from sys import argv
import networkx as nx
import pylab
import numpy as np
from itertools import combinations
from steiner_midpoint import slope_vector, delta_point
from numpy.random import exponential
from dist_functions import *

SYNAPSE_RATE = 0.1

def new_synapse_points(coord1, coord2, rate=SYNAPSE_RATE, offset=0):
    slope_vec = slope_vector(coord1, coord2)
    #max_dist = point_dist(coord1, coord2)
    length = point_dist(coord1, coord2)
    npoints = int((length - offset) / rate)
    new_offset = (length - offset) - (npoints * rate)
    assert 0 <= new_offset <= rate
    new_offset = rate - new_offset
    assert 0 <= new_offset <= rate
    dist = 0
    beta = 1.0 / rate
    new_points = []
    d = offset / length
    for i in xrange(npoints):
        new_point = delta_point(coord1, slope_vec, d)
        new_points.append(new_point)
        d += rate / length
    return new_points

def add_synapses(G, rate=SYNAPSE_RATE):
    G.graph['synapses'] = []
    next_node = max(G.nodes()) + 1
    for u, v in G.edges():
        new_points = new_synapse_points(G.node[u]['coord'], G.node[v]['coord'],\
                                        rate=rate)
        new_nodes = []
        for new_point in new_points:
            new_nodes.append(next_node)
            G.add_node(next_node)
            G.node[next_node]['coord'] = new_point
            G.node[next_node]['label'] = 'synapse'
            G.graph['synapses'].append(next_node)
            next_node += 1

        segment_nodes = [u] + new_nodes + [v]
        G.remove_edge(u, v)
        for i in xrange(len(segment_nodes) - 1):
            n1 = segment_nodes[i]
            n2 = segment_nodes[i + 1]
            G.add_edge(n1, n2)
            G[n1][n2]['length'] = node_dist(G, n1, n2)
            print G[n1][n2]['length']

def get_neuron_points(filename, dim='3D'):
    #for arbor_type in ["2","3","4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
    graphs = []
    for arbor_type in ["2", "3", "4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
        
        # Reads in 3D arborization.
        G = nx.Graph()
        P2Coord = {} # u-> (x,y)
        root = -1
        with open(filename) as f:
            for line in f:
                if line.startswith("#"): continue

                cols = line.strip().split()
                assert len(cols) == 7

                if not (cols[1] == "1" or cols[1] == arbor_type): continue
 
                if cols[6] == "-1":
                    if root != -1: assert False # duplicate root.

                    root = int(cols[0])
                                        
                    assert root not in P2Coord and root not in G
                    G.add_node(root)
                    coord = None
                    if dim == "3D":
                        coord = (float(cols[2]),float(cols[3]),float(cols[4]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        coord = (float(cols[2]),float(cols[3]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]))
                    else:
                        assert False
                    #G.node[root]['coord'] = pylab.array(coord)
                    G.node[root]['coord'] = coord
                    G.graph['root'] = root


                else:
                    u,v = int(cols[6]),int(cols[0])
                    assert u in G #and u in P2Coord
                    assert not G.has_edge(u,v)
            
                    coord = None
                    if dim == "3D":
                        coord = (float(cols[2]),float(cols[3]),float(cols[4]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        coord = (float(cols[2]),float(cols[3]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]))
                    else:
                        assert False
                    G.add_edge(u, v)
                    
                    #G.node[v]['coord'] = pylab.array(coord)
                    G.node[v]['coord'] = coord
                    
                    G[u][v]['length'] = point_dist(G.node[u]['coord'], G.node[v]['coord'])
        
        if G.number_of_nodes() > 0:
            assert 'root' in G.graph
            label_points(G)
            #add_synapses(G)
            graphs.append(G)
        else:
            graphs.append(None)

    return graphs

def initialize_lengths(G):
    for u, v in G.edges_iter():
        p1, p2 = G.node[u]['coord'], G.node[v]['coord']
        G[u][v]['length'] = point_dist(p1, p2)

def get_label(G, u):
    if 'label' in G.node[u] and 'label' not in ['root', 'tip', 'branch', 'continue']:
        return G.node[u]['label']
    if u == G.graph['root']: 
        return "root"
    elif G.degree(u) == 1: 
        return "tip"
    elif G.degree(u) > 2: 
        return "branch"
    else:                  
        return "continue"

def label_points(G):
    root = G.graph['root']
    for u in G:
        G.node[u]['label'] = get_label(G, u)

def viz_tree(G, name, outdir='figs'):
    """ Displays plant/tree visualization. """
    
    root = G.graph['root']
    node_size,node_color = [],[]
    pos = {}    
    label_points(G)
    for u in G:
        coord = G.node[u]['coord']
        pos[u] = (coord[0], coord[1])
        
        label = None
        if 'label' in G.node[u]:
            label = G.node[u]['label']
        else:
            label = get_label(G, u)
        if label == "root":
            node_color.append('black')
            node_size.append(350)
        elif label == 'synapse':
            node_color.append('green')
            node_size.append(15)
        elif label == 'tip':
            node_color.append('blue')
            node_size.append(20)
        elif label == 'branch':
            node_color.append('blue')
            node_size.append(10)
        elif label == "continue":
            node_color.append('blue')
            #node_size.append(250)   
            node_size.append(10)   
        elif label == 'centroid':
            node_color.append('purple')
            node_size.append(250)
        elif label == 'isolated_start':
            node_color.append('yellow')
            node_size.append(100)
        elif label == 'isolated_end':
            node_color.append('blue')
            node_size.append(100)
        elif label == 'steiner_midpoint':
            node_color.append('blue')
            node_size.append(1)
        else:
            print label
            assert False


    nx.draw(G,pos=pos,arrows=False,with_labels=False,node_size=node_size,node_color=node_color,edge_color="brown",width=4,font_size=12,font_color='red',font_weight='bold')
    #nx.draw(G,pos=pos,arrows=False,with_labels=False,node_size=15,font_size=12,font_color='red',font_weight='bold')
    pylab.draw()
    #PP.show()
    
    '''
    mcost = mst_cost(G)
    scost = satellite_cost(G)
    title_str = 'satellite cost = %f\nspanning tree cost = %f' % (scost, mcost)
    #pylab.text(125, 125, title_str)
    '''

    pylab.savefig("%s/%s.pdf" % (outdir, name))
    pylab.close()


def main():
    #neuron_file = argv[1]
    #print map(str, neuron_info(neuron_file))
    pass

if __name__ == '__main__':
    main()
