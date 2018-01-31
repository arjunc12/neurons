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
from graph_utils import is_tree

SYNAPSE_RATE = 1.8661051308416368

def truncate_graph(G):
    remove_nodes = set()
    queue = [G.graph['root']]
    new_root = None
    while len(queue) > 0:
        assert len(queue) == 1
        curr = queue.pop()
        degree = 0
        descendants = []
        for u in G.neighbors(curr):
            if u not in remove_nodes:
                degree += 1
                descendants.append(u)
        if degree <= 1:
            queue += descendants
            remove_nodes.add(curr)
        else:
            done = True
            new_root = curr

    G.remove_nodes_from(remove_nodes)
    G.graph['root'] = new_root

def new_synapse_points(coord1, coord2, rate=SYNAPSE_RATE, offset=0):
    #assert 0 <= offset <= rate
    slope_vec = slope_vector(coord1, coord2)
    max_dist = point_dist(coord1, coord2)
    total_dist = offset
    dists = []
    while total_dist < max_dist:
        dists.append(total_dist)
        dist = exponential(rate)
        total_dist += dist
    new_offset = total_dist - max_dist
    new_points = []
    for dist in dists:
        assert 0 <= dist <= max_dist
        d = dist / max_dist
        new_point = delta_point(coord1, slope_vec, d)
        new_points.append(new_point)
    return new_points, new_offset

def add_synapses(G, rate=SYNAPSE_RATE):
    H = G.copy()
    assert is_tree(H)
    root = H.graph['root']
    stack = [root]
    offsets = {}
    offsets[root] = 0
    visited = set()
    H.graph['synapses'] = []
    next_node = max(G.nodes()) + 1
    while len(stack) > 0:
        curr = stack.pop()
        visited.add(curr)
        for n in H.neighbors(curr):
            if n in visited:
                continue

            coord1 = H.node[curr]['coord']
            coord2 = H.node[n]['coord']
            offset = offsets[curr]
            new_points, new_offset = new_synapse_points(coord1, coord2,\
                                                        rate=rate,\
                                                        offset=offset)
            
            new_nodes = []
            for new_point in new_points:
                new_nodes.append(next_node)
                H.add_node(next_node)
                H.node[next_node]['coord'] = new_point
                H.node[next_node]['label'] = 'synapse'
                H.graph['synapses'].append(next_node)
                next_node += 1

            segment_nodes = [curr] + new_nodes + [n]
            H.remove_edge(curr, n)
            for i in xrange(len(segment_nodes) - 1):
                n1 = segment_nodes[i]
                n2 = segment_nodes[i + 1]
                H.add_edge(n1, n2)
                H[n1][n2]['length'] = node_dist(H, n1, n2)

            stack.append(n)
            offsets[n] = new_offset

    print len(H.graph['synapses']), 'synapses'
    return H

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

    if graphs[0] != None:
        G = graphs[0].copy()
        truncate_graph(G)
        if G.number_of_nodes() > 0:
            assert 'root' in G.graph
            label_points(G)
            graphs.append(G)
        else:
            graphs.append(None)
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
            node_size.append(10)
        elif label == 'branch':
            node_color.append('blue')
            node_size.append(10)
        elif label == "continue":
            node_color.append('brown')
            #node_size.append(250)   
            node_size.append(1)
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
            node_color.append('brown')
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
    G = nx.Graph()
    for i in xrange(1, 11):
        G.add_edge(i - 1, i)
    print G.nodes()
    print G.edges()
    for j in xrange(11, 20):
        pass
        #G.add_edge(10, j)
    G.graph['root'] = 0
    truncate_graph(G)
    print G.nodes()
    print G.edges()

if __name__ == '__main__':
    main()
