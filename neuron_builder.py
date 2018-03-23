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

from random import uniform, random, expovariate, gauss, choice, seed

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
    while G.degree(curr) == 1 and G.node[curr]['label'] not in ['synapse', 'root']:
        parent = G.neighbors(curr)[0]
        G.remove_node(curr)
        curr = parent
            
def retract_graph(G):
    leaves = []
    for u in G.nodes_iter():
        if (G.node[u]['label'] not in ['synapse', 'root']) and (G.degree(u) == 1):
            leaves.append(u)
    
    for leaf in leaves:
        retract_branch(G, leaf)
        
def normalize_vector(vector):
    total = 0
    for val in vector:
        total += val ** 2
    total **= 0.5
    norm_vec = []
    for val in vector:
        norm_vec.append(val / total)
    return tuple(norm_vec)
        
def extend(coord1, coord2, trial_length=None):
    slope = normalize_vector(slope_vector(coord1, coord2))
    if trial_length == None:
        trial_length = expovariate(1)
    return delta_point(coord2, slope, trial_length)
    
def add_extension_node(G, u, parent, trial_length=None):
    coord1 = G.node[parent]['coord']
    coord2 = G.node[u]['coord']
    new_coord = extend(coord1, coord2, trial_length=trial_length)
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
        
def add_bifurcations(G, u, dim=3, bifurcations=2, dist=None):
    coord = G.node[u]['coord']
    new_nodes = []
    if dist == None:
        dist = expovariate(1)
    for i in xrange(bifurcations):
        next_coord = None
        done = False
        while not done:
            next_coord = []
            magnitude = 0
            for i in xrange(dim):
                delta = gauss(0, 1)
                next_coord.append(delta)
                magnitude += delta ** 2
            magnitude **= 0.5
            if magnitude == 0:
                continue
            for i in xrange(len(next_coord)):
                delta = next_coord[i]
                next_coord[i] = coord[i] + (dist * delta / magnitude)
            pdist = point_dist(coord, next_coord)
            if pdist == dist:
                done = True
            
        next_coord = tuple(next_coord)
        new_nodes.append(add_new_node(G, next_coord, u))
    return new_nodes

def annihilate(G, u, detection_radius):
    coord1 = G.node[u]['coord']
    for v in G.nodes_iter():
        if v == u:
            continue
        coord2 = G.node[v]['coord']
        if point_dist(coord1, coord2) <= detection_radius:
            return True
    return False

def can_extend(G, unmarked_points, trial_length, r_puncta):
    for u in G.nodes_iter():
        coord1 = G.node[u]['coord']
        for coord2 in unmarked_points:
            if point_dist(coord1, coord2) <= trial_length + r_puncta:
                return True
    return False 

def update_graph_barw(G, unmarked_points, dim=3, **kwargs):
    branching_prob = kwargs['branching_prob']
    trial_length = kwargs['trial_length']
    r_annihilation = kwargs['radius_annihilation']
    r_puncta = kwargs['radius_puncta']
    for u in G.nodes():
        if u == G.graph['root']:
            if random() <= branching_prob:
                add_bifurcations(G, u, dim=dim, bifurcations=1)
        elif G.degree(u) == 1 and G.node[u]['state'] == 'active':
            if annihilate(G, u, r_annihilation):
                G.node[u]['state'] = 'stabilized'
            connect_to_synapses(G, u, unmarked_points, radius=1)
            bifurcate = random() <= branching_prob
            if bifurcate:
                add_bifurcations(G, u, dim=dim, bifurcations=1)
            else:
                parent = G.neighbors(u)[0]
                add_extension_node(G, u, parent, trial_length)
                
    return True
                
def update_graph_snider(G, unmarked_points, dim=3, **kwargs):
    r_puncta = kwargs['radius_puncta']
    trial_length = kwargs['trial_length']
    r_remove = kwargs['radius_remove']
    u = choice(G.nodes())
    new_extension = add_bifurcations(G, u, dim=dim, bifurcations=1, dist=trial_length)[0]
    new_synapses = connect_to_synapses(G, new_extension, unmarked_points, radius=r_puncta)
    if len(new_synapses) == 0:
        G.remove_node(new_extension)
        
    return can_extend(G, unmarked_points, trial_length, r_puncta)

def get_update_func(algorithm):
    if algorithm == 'snider':
        return update_graph_snider
    elif algorithm == 'barw':
        return update_graph_barw
    else:
        raise ValueError("Invalid Algorithm")

def init_graph(dim=3):
    G = nx.Graph()
    G.graph['synapses'] = []
    G.add_node(1)
    G.node[1]['coord'] = tuple([0] * dim)
    G.graph['root'] = 1
    G.node[1]['label'] = 'root'
    return G

def update_graph(G, algorithm, unmarked_points, dim=3, **kwargs):
    update_func = get_update_func(algorithm)
    return update_func(G, unmarked_points, dim=dim, **kwargs)
                
def build_neuron(algorithm='snider', dim=3, **kwargs):
    unmarked_points = grid_points3d(xmin=-2, xmax=2, ymin=-2, ymax=2, zmin=0, zmax=0)
    G = init_graph(dim=dim)
    done = False
    while not done:
        if len(unmarked_points) == 0:
            break
        can_extend = update_graph(G, algorithm, unmarked_points, dim=dim, **kwargs)
        done = not can_extend
        
    retract_graph(G)
    print G.number_of_nodes(), G.number_of_edges()

def build_neuron_video(algorithm='snider', dim=3, **kwargs):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    init_points = grid_points3d(xmin=-2, xmax=2, ymin=-2, ymax=2, zmin=0, zmax=0)
    points = init_points[:]
    G = init_graph(dim=dim)
    graphs = [G.copy()]
    unmarked_points = [init_points]

    done = False
    while not done:
        if len(points) == 0:
            break
        can_extend = update_graph(G, algorithm, points, dim=dim, **kwargs)
        graphs.append(G.copy())
        unmarked_points.append(points[:])
        done = not can_extend
        
    retract_graph(G)
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
    parser.add_argument('-rp', '--radius_puncta', default=1)
    parser.add_argument('-rr', '--radius_remove', default=1)
    parser.add_argument('-ra', '--radius_annihilation', default=1)
    parser.add_argument('-p', '--branching_prob', default=0.1)
    parser.add_argument('-l', '--trial_length', default=1)
    parser.add_argument('-v', '--video', action='store_true')
    
    args = parser.parse_args()
    algorithm = args.algorithm
    dim = args.dim
    video = args.video
    kwargs = vars(args)
    del kwargs['algorithm'], kwargs['dim']

    if video:
        build_neuron_video(algorithm, dim, **kwargs)
    else:
        build_neuron(algorithm, dim, **kwargs)
    #G = init_graph()
    #add_bifurcations(G, 1, dim=3, bifurcations=1, dist=0.775567698885)

if __name__ == '__main__':
    main()
