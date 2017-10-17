import matplotlib as mpl
mpl.use('agg')
from sys import argv
import networkx as nx
import pylab
import numpy as np
from itertools import combinations
from scipy.spatial.distance import euclidean

def pareto_cost(mcost, scost, alpha):
    # alpha = 0 minimizes satellite cost
    # alpha = 1 minimizes spanning tree cost
    mcost *= alpha
    scost *= (1 - alpha)
    cost = mcost + scost
    return cost 

def centroid(G):
    root = G.graph['root']
    root_coord = G.node[root]['coord']
    centroid = np.zeros(len(root_coord))
    for u in G.nodes_iter():
        point = G.node[u]['coord']
        assert len(point) == len(root_coord)
        if u != root:
            centroid += point
    centroid /= G.number_of_nodes() - 1
    return centroid

def new_synapse_points(coord1, coord2, dist=1):
    new_points = []
    pass
    return new_points

def add_synapses(G, dist=1):
    latest_node = max(G.nodes())
    for u, v in G.edges():
        new_points = add_synapse_points(G.node[u]['coord'], G.node[v]['coord'])
        new_nodes = []
        for new_point in new_point:
            new_node = latest_node + 1
            new_nodes.append(new_node)
            G.add_node(new_node)
            G.node[new_node]['coord'] = new_point
            latest_node = new_node

        new_nodes = [u] + new_nodes + [y]
        for i in xrange(len(new_nodes) - 1):
            x, y = new_nodes[i], new_nodes[i + 1]
            pass

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

        assert 'root' in G.graph
        label_points(G)
        graphs.append(G)

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
        elif label in ['branch', 'tip']:
            node_color.append('green')
            node_size.append(20)
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
        elif label == 'steiner_point':
            node_color.append('blue')
            node_size.append(20)
        else:
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

def point_dist(p1, p2):
    '''
    '''
    assert len(p1) == len(p2)
    sq_dist = 0
    for i in xrange(len(p1)):
        x1, x2 = p1[i], p2[i]
        sq_dist += (x1 - x2) ** 2
    return sq_dist ** 0.5
    '''
    '''
    #return euclidean(p1, p2)
    #return np.linalg.norm(p1 - p2)

def node_dist(G, u, v):
    p1 = G.node[u]['coord']
    p2 = G.node[v]['coord']
    return point_dist(p1, p2)

def root_dist(G, u):
    assert 'root' in G.graph
    return node_dist(G, u, G.graph['root'])

def sort_neighbors(G):
    for u in G.nodes_iter():
        G.node[u]['close_neighbors'] = sorted(G.neighbors(u), key = lambda v : G[u][v]['length'])
    G.graph['sorted'] = True

def pareto_dist(pareto_mcosts, pareto_scosts, mcost, scost):
    best_dist = float("inf")
    best_index = None

    assert len(pareto_mcosts) == len(pareto_scosts)

    for i in xrange(len(pareto_mcosts)):
        pareto_mcost = pareto_mcosts[i]
        pareto_scost = pareto_scosts[i]

        dist = point_dist((pareto_mcost, pareto_scost), (mcost, scost))
        if dist < best_dist:
            best_dist = dist
            best_index = i

    return best_dist, best_index 

def main():
    #neuron_file = argv[1]
    #print map(str, neuron_info(neuron_file))
    pass

if __name__ == '__main__':
    main()
