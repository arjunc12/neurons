import os
from neuron_utils import get_neuron_points, viz_tree
from graph_utils import is_tree
from sys import argv

import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from random_graphs import random_mst, random_point_graph, random_points

from steiner_midpoint import slope_vector, delta_point

import networkx as nx

from random import uniform, random

from dist_functions import point_dist

BUILDER_EXE = 'neuronbuilder2'

OUTDIR = 'neuron_builder'
OUTFILE = 'neuron_builder.swc'

def build_neuron_snider(radius=1):
    build_command = './%s %f > %s' % (BUILDER_EXE, radius, OUTFILE)
    print build_command
    os.system(build_command)
    G = get_neuron_points(OUTFILE)[1]
    assert is_tree(G)
    return G

def draw_graph_3d(G, ax):
    for u in G.nodes_iter():
        coord = G.node[u]['coord']
        x = [coord[0]]
        y = [coord[1]]
        z = [coord[2]]
        ax.scatter(x, y, z, c='r', s=100)
        
    for u, v in G.edges_iter():
        coord1, coord2 = G.node[u]['coord'], G.node[v]['coord']
        x = [coord1[0], coord2[0]]
        y = [coord1[1], coord2[1]]
        z = [coord1[2], coord2[2]]
        ax.plot(x, y, z, c='k')

def add_new_node(G, new_coord, u, synapse=False):
    next_node = G.number_of_nodes() + 1
    G.add_node(next_node)
    G.node[next_node]['coord'] = new_coord
    G.add_edge(u, next_node)
    if synapse:
        G.graph['synapses'].append(next_node)
        G.node[next_node]['label'] = 'synapse'
    else:
        G.node[next_node]['label'] = 'steiner_midpoint'
        
def next_choice():
    choice = random()
    if choice <= 0.2:
        return 'extend'
    elif 0.2 < choice <= 0.25:
        return "bifurcate" 
        
def extend(coord1, coord2):
    slope = slope_vector(coord1, coord2)
    dist = uniform(1, 2)
    return delta_point(coord2, slope, dist)
    
def add_extension_node(G, u, parent):
    coord1 = G.node[parent]['coord']
    coord2 = G.node[u]['coord']
    new_coord = extend(coord1, coord2)
    add_new_node(G, new_coord, u)

def connect_to_synapses(G, u, unmarked_points, radius=1):
    coord = G.node[u]['coord']
    added_points = []
    for point in unmarked_points:
        if point_dist(coord, point) <= radius:
            added_points.append(point)
    
    for point in added_points:
        add_new_node(G, point, u, synapse=True)
        unmarked_points.remove(point)

def add_bifurcations(G, u, n=2):
    coord = G.node[u]['coord']
    for i in xrange(n):
        next_coord = []
        for i in xrange(3):
            next_coord.append(uniform(-1, 1))
        next_coord = tuple(next_coord)
        add_new_node(G, next_coord, u)       

def update_graph(G, unmarked_points, radius=1):
    for u in G.nodes():
        if u == G.graph['root']:
            add_bifurcations(G, u, n=1)
        elif G.degree(u) == 1:
            connect_to_synapses(G, u, marked_points, unmarked_points, radius=1)
            choice = next_choice()
            if choice == 'extend':
                parent = G.neighbors(u)[0]
                add_extension_node(G, u, parent)
            elif choice == 'bifurcate':
                add_bifurcations(G, u)
                
def build_neuron_video():
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    unmarked_points = set(random_points())
    G = nx.Graph()
    G.graph['synapses'] = []
    G.add_node(1)
    G.node[1]['coord'] = (0, 0, 0)
    G.graph['root'] = 1
    prev_node = None
    graphs = [G.copy()]
    for i in xrange(100):
        update_graph(G, unmarked_points, radius=1)
        graphs.append(G.copy())
        

    def init():
        x, y, z = zip(*unmarked_points)
        ax.scatter(x, y, zs=z, s=250)
        
    def redraw(frame):
        ax.clear()
        init()
        draw_graph_3d(graphs[frame], ax)
        
    def redraw2d(frame):
        x, y, z = zip(*unmarked_points)
        plt.clf()
        plt.scatter(x, y, s=250)
        pos = {}
        G = graphs[frame]
        for u in G.nodes_iter():
            coord = G.node[u]['coord']
            pos[u] = (coord[0], coord[1])
        nx.draw_networkx(G, pos=pos, with_labels=False)
        plt.draw()
        
    ani = animation.FuncAnimation(fig, redraw2d, init_func=init, frames=len(graphs), \
                                  interval = 1000)
    ani.save('neuron_builder/neuron_builder.mp4')
            
    plt.close()
        

def main():
    build_neuron_video()

if __name__ == '__main__':
    main()
