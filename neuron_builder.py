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

from random import uniform, random, expovariate, gauss, choice

from dist_functions import point_dist

import argparse

BUILDER_EXE = 'neuronbuilder2'

OUTDIR = 'neuron_builder'
OUTFILE = 'neuron_builder.swc'

def grid_points2d(xmin=-10, xmax=10, ymin=-10, ymax=10, dist=1):
    x = xmin
    y = ymin
    points = []
    
    i = 0
    while x <= xmax:
        while y <= ymax:
            points.append((x, y))
            y += dist
        x += dist
    return points
    
def grid_points3d(xmin=-10, xmax=10, ymin=-10, ymax=10, zmin=-10, zmax=10, dist=1):
    points = []    
    x = xmin
    while x <= xmax:
        y = ymin
        while y <= ymax:
            z = zmin
            while z <= zmax: 
                points.append((x, y, z))
                z += dist
            y += dist
        x += dist
    return points               

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
    G.node[next_node]['state'] = 'active'
    return next_node
    
def retract_branch(G, u):
    curr = u
    while G.degree(curr) == 1:
        parent = G.neighbors(curr)[0]
        G = nx.contracted_edge(G, (parent, curr))
        curr = parent
    return G
    
def retract_graph(G):
    leaves = []
    for u in G.nodes_iter():
        if (G.node[u]['label'] != 'synapse') and (G.degree(u) == 1):
            leaves.append(u)
    
    print G.number_of_nodes()
    for leaf in leaves:
        G = retract_branch(G, leaf)
        
    return G
            
def next_choice():
    choice = random()
    if choice <= 0.2:
        return 'stabilize'
    elif 0.2 < choice <= 0.25:
        return "bifurcate"
    else:
        return "extend"

def normalize_slope(slope_vec):
    total = 0
    for val in slope_vec:
        total += val ** 2
    total **= 0.5
    norm_slope = []
    for val in slope_vec:
        norm_slope.append(val / total)
    return tuple(norm_slope)
        
def extend(coord1, coord2):
    slope = normalize_slope(slope_vector(coord1, coord2))
    dist = expovariate(1)
    return delta_point(coord2, slope, dist)
    
def add_extension_node(G, u, parent):
    coord1 = G.node[parent]['coord']
    coord2 = G.node[u]['coord']
    new_coord = extend(coord1, coord2)
    return add_new_node(G, new_coord, u)

def connect_to_synapses(G, u, unmarked_points, radius=1):
    coord = G.node[u]['coord']
    added_points = []
    for point in unmarked_points:
        if point_dist(coord, point) <= radius:
            added_points.append(point)

    new_nodes = []    
    for point in added_points:
        new_nodes.append(add_new_node(G, point, u, synapse=True))
        unmarked_points.remove(point)
    return new_nodes
        
def add_bifurcations(G, u, dim=3, bifurcations=2):
    coord = G.node[u]['coord']
    new_nodes = []
    for i in xrange(bifurcations):
        next_coord = []
        for i in xrange(dim):
            next_coord.append(coord[i] + gauss(0, 0.1))
        next_coord = tuple(next_coord)
        new_nodes.append(add_new_node(G, next_coord, u))
    return new_nodes

def update_graph_barw(G, unmarked_points, dim=3, **kwargs):
    detection_radius = kwargs['detection_radius']
    branching_prob = kwargs['branching_prob']
    trial_length = kwargs['trial_length']
    for u in G.nodes():
        if u == G.graph['root']:
            if random() <= 0.25:
                add_bifurcations(G, u, dim=dim, bifurcations=1)
        elif G.degree(u) == 1 and G.node[u]['state'] == 'active':
            connect_to_synapses(G, u, unmarked_points, radius=1)
            choice = next_choice()
            if choice == 'extend':
                parent = G.neighbors(u)[0]
                add_extension_node(G, u, parent)
            elif choice == 'bifurcate':
                add_bifurcations(G, u, dim=dim, bifurcations=1)
            elif choice == 'stabilize':
                G.node[u]['state'] = 'stabilized'
                
def update_graph_snider(G, unmarked_points, dim=3, **kwargs):
    detection_radius = kwargs['detection_radius']
    trial_length = kwargs['trial_length']
    u = choice(G.nodes())
    new_extension = add_bifurcations(G, u, dim=dim, bifurcations=1)[0]
    new_synapses = connect_to_synapses(G, new_extension, unmarked_points, radius=radius)
    if len(new_synapses) == 0:
        G.remove_node(new_extension)

def get_update_func(algorithm):
    if algorithm == 'snider':
        return update_graph_sinder
    elif algorithm == 'barw':
        return update_graph_barw
    else:
        raise ValueError("Invalid Algorithm")

def update_graph(G, algorithm, unmarked_points, dim=3, **kwargs):
    update_func = get_update_func(algorithm)
    update_graph(G, unmarked_points, dim=dim, **kwargs)
                
def build_neuron_video(algorithm='snider', dim=3, **kwargs):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    init_points = grid_points3d(xmin=-2, xmax=2, ymin=-2, ymax=2, zmin=0, zmax=0)
    points = init_points[:]
    
    G = nx.Graph()
    G.graph['synapses'] = []
    G.add_node(1)
    G.node[1]['coord'] = tuple([0] * dim)
    G.graph['root'] = 1
    G.node[1]['label'] = 'root'
    graphs = [G.copy()]
    unmarked_points = [init_points]

    for i in xrange(100):
        if len(points) == 0:
            break
        
        graphs.append(G.copy())
        unmarked_points.append(points[:])
        
    G = retract_graph(G)
    graphs.append(G.copy())
    unmarked_points.append(points[:])

    def init():
        x, y, z = zip(*init_points)
        ax.scatter(x, y, zs=z, s=25)
        
    def redraw(frame):
        plt.clf()
        G = graphs[frame]
        for u in G.nodes_iter():
            if G.node[u]['label'] == 'synapse':
                points.remove(G.node[u]['coord'])
        init()            
        draw_graph_3d(graphs[frame], ax)
        
    def redraw2d(frame):
        print frame
        plt.clf()
        G = graphs[frame]
        if len(unmarked_points[frame]) > 0:
            x, y, z = zip(*unmarked_points[frame])
            plt.scatter(x, y, s=250)
        pos = {}
        for u in G.nodes_iter():
            coord = G.node[u]['coord']
            pos[u] = (coord[0], coord[1])
        viz_tree(G, save=False)
        
    ani = animation.FuncAnimation(fig, redraw2d, init_func=init, frames=len(graphs), \
                                  interval = 1000)
    ani.save('neuron_builder/neuron_builder.mp4')
            
    plt.close()
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--algorithm', default='snider')
    parser.add_argument('-d', '--dim', default=3)
        
    
    args = parser.parse_args()
    algorithm = args.algorithm
    dim = args.dim
    kwargs = vars(args)
    del kwargs['algorithm']
    del kwargs['dim']

    build_neuron_video()

if __name__ == '__main__':
    main()
