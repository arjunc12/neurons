import os
from neuron_utils import get_neuron_points, viz_tree, get_neuron_points, retract_graph
from steiner_midpoint import slope_vector, best_midpoint_approx
from graph_utils import is_tree
from sys import argv

import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from random_graphs import *

from steiner_midpoint import slope_vector, delta_point, steiner_points

import networkx as nx

from random import uniform, random, expovariate, gauss, choice, seed

from dist_functions import point_dist, node_dist, line_seg_dist

import math

import argparse

MAX_STEPS = 1000

BUILDER_EXE = 'neuronbuilder2'

OUTDIR = '/iblsn/data/Arjun/neurons/neuron_builder'
OUTFILE = 'neuron_builder.swc'

STEINER_MIDPOINTS = 10

MIN_ANGLE = -math.pi / 3
MAX_ANGLE = math.pi / 3

#seed(10)

def grid_points2d(xmin=-10, xmax=10, ymin=-10, ymax=10, dist=1):
    x = xmin
    y = ymin
    points = []
    
    i = 0
    while x <= xmax:
        while y <= ymax:
            if (x != 0) or (y != 0):
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
                if (x != 0) or (y != 0) or (z != 0):
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
    for u in G.nodes():
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

def remove_close_points(point, unmarked_points, remove_radius=1):
    close_points = []
    for unmarked_point in unmarked_points:
        if point_dist(point, unmarked_point) <= remove_radius:
            close_points.append(unmarked_point)
    for close_point in close_points:
        unmarked_points.remove(close_point)

def connect_to_synapses(G, u, unmarked_points, puncta_radius=1, remove_radius=1):
    coord = G.node[u]['coord']
    added_points = []
    for point in unmarked_points:
        if point_dist(coord, point) <= puncta_radius:
            added_points.append(point)
    
    if len(added_points) == 0:
        return None
    
    new_point = choice(added_points)
    remove_close_points(new_point, unmarked_points, remove_radius=remove_radius)
    new_node = add_new_node(G, new_point, u, synapse=True) 
    return new_node

def self_crossing(G, P0, Q1, radius_annihilation):
    for Q0, Q1 in G.graph['line segments']:
        dist = line_seg_dist(P0, P1, Q0, Q1)
        if dist <= radius_annihilation:
            return True
    return False
 
def normalize_vector(vector):
    magnitude = 0
    for delta in vector:
        magnitude += delta ** 2
    magnitude **= 0.5
    for i in xrange(len(vector)):
        vector[i] /= magnitude
        
def add_bifurcations(G, u, dim=3, bifurcations=2, dist=None, **kwargs):
    coord_mins = None
    coord_maxes = None
    if 'coord_mins' in kwargs:
        coord_mins = kwargs['coord_mins']
    if 'coord_maxes' in kwargs:
        coord_maxes = kwargs['coord_maxes']
    
    coord = G.node[u]['coord']
    assert len(coord) == dim
    
    parent = G.node[u]['parent']
    parent_coord = None
    if parent != None:
        parent_coord = G.node[parent]['coord']
    else:
        parent_coord = []
        for i in xrange(dim):
            parent_coord.append(gauss(0, 1))
    assert len(parent_coord) == dim
    parent_coord = tuple(parent_coord)
      
    slope = slope_vector(parent_coord, coord)
    normalize_vector(slope)
    slope_variance = 0.1
        
    new_nodes = []
    for i in xrange(bifurcations):
        next_coord = []
        
        if dist == None:
            dist = expovariate(1)
        
        for i in xrange(dim):
            delta = slope[i] + gauss(0, slope_variance)
            next_coord.append(coord[i] + (dist * delta))
        
        for i in xrange(len(next_coord)):
            if coord_mins != None:
                next_coord[i] = max(next_coord[i], coord_mins[i])
            if coord_maxes != None:
                next_coord[i] = min(next_coord[i], coord_maxes[i])
            
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
    for u in G.nodes():
        coord1 = G.node[u]['coord']
        for coord2 in unmarked_points:
            if point_dist(coord1, coord2) <= trial_length + r_puncta:
                return True
    return False 

def update_graph_barw(G, unmarked_points, alpha, dim=3, **kwargs):
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

def add_midpoints(G, u, v):
    G.remove_edge(u, v)
    node_index = G.number_of_nodes() + 1
    
    p1 = G.node[u]['coord']
    p2 = G.node[v]['coord']
   
    midpoints = steiner_points(p1, p2, npoints=STEINER_MIDPOINTS)
    midpoint_nodes = []
    for midpoint in midpoints:
        midpoint_node = node_index
        G.add_node(midpoint_node)
        G.node[midpoint_node]['coord'] = midpoint
        midpoint_nodes.append(midpoint_node)
        node_index += 1

    line_nodes = [v] + list(reversed(midpoint_nodes)) + [u]
    for i in xrange(-1, -len(line_nodes), -1):
        n1 = line_nodes[i]
        n2 = line_nodes[i - 1]
        G.add_edge(n1, n2)
        G[n1][n2]['length'] = node_dist(G, n1, n2)
        if 'label' not in G.node[n2]:
            G.node[n2]['label'] = 'steiner_midpoint'
        G.node[n2]['parent'] = n1
    G.graph['line segments'].append((u, v))

def check_parents(G):
    for u in G.nodes():
        assert 'parent' in G.node[u]
        parent = G.node[u]['parent']
        if u == G.graph['root']:
            assert parent == None
        else:
            assert G.has_edge(u, parent)
                
def update_graph_snider(G, unmarked_points, alpha, dim=3, **kwargs):
    trial_length = kwargs['trial_length']
    r_puncta = kwargs['radius_puncta']
    r_remove = kwargs['radius_remove']
    p_root = kwargs['prob_root']
    
    coord_mins = [kwargs['xmin'], kwargs['ymin'], kwargs['zmin']]
    coord_maxes = [kwargs['xmax'], kwargs['ymax'], kwargs['zmax']]
    
    u = None
    root = G.graph['root']
    candidates = list(G.nodes())
    assert root in candidates
    if random() <= p_root or len(candidates) == 1:
        u = G.graph['root']
    else:
        candidates.remove(root)
        u = choice(candidates)
    new_extension = add_bifurcations(G, u, dim=dim, bifurcations=1, dist=trial_length,\
                                     coord_mins=coord_mins, coord_maxes=coord_maxes)
    new_extension = new_extension[0]
    new_synapse = connect_to_synapses(G, new_extension, unmarked_points,\
                                      puncta_radius=r_puncta, remove_radius=r_remove)
    
    if new_synapse == None:
        G.remove_node(new_extension)
    else:
        p1 = G.node[u]['coord']
        p2 = G.node[new_extension]['coord']
        p3 = G.node[new_synapse]['coord']
        best_midpoint, best_choice = best_midpoint_approx(p1, p2, p3, alpha)
        if best_choice == 1:
            G.remove_node(new_extension)
            G.add_edge(u, new_synapse)
            nx.relabel_nodes(G, {new_synapse : new_extension}, copy=False)
            new_synapse = new_extension
            new_extension = u
        else:
            G.node[new_extension]['coord'] = best_midpoint
            add_midpoints(G, u, new_extension)
        add_midpoints(G, new_extension, new_synapse)
        
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
    G.node[1]['parent'] = None
    G.graph['line segments'] = []
    return G

def update_graph(G, algorithm, unmarked_points, alpha, dim=3, **kwargs):
    update_func = get_update_func(algorithm)
    return update_func(G, unmarked_points, alpha, dim=dim, **kwargs)

def swc_line(G, u, parents, point_labels):
    write_items = []
    root = u == G.graph['root']
    
    next_label = u
    point_labels[u] = next_label
    write_items.append(next_label)
    
    segment_type = 3
    write_items.append(segment_type)

    write_items += list(G.node[u]['coord'])
    write_items.append(0)
    if root:
        write_items.append(-1)
    else:
        assert u in parents
        parent = parents[u]
        assert parent in point_labels
        write_items.append(point_labels[parent])
    write_items = map(str, write_items)
    write_items = ' '.join(write_items)
    return write_items

def write_to_swc(G, outfile='neuron_builder.swc'):
    with open(outfile, 'w') as f:
        root = G.graph['root']
        queue = [root]
        visited = set()
        point_labels = {}
        parents = {}
        while len(queue) > 0:
            curr = queue.pop(0)
            line = swc_line(G, curr, parents, point_labels)
            f.write('%s\n' % line)
            visited.add(curr)
            for n in G.neighbors(curr):
                if n not in visited:
                    parents[n] = curr
                    queue.append(n)

def read_tree(tree_dir):
    graphs = get_neuron_points('%s/tree.swc' % tree_dir)
    G = graphs[1]
    G.graph['synapses'] = []
    for line in open('%s/synapses.txt' % tree_dir):
        synapse = int(line)
        G.node[synapse]['label'] = 'synapse'
        G.graph['synapses'].append(synapse)
    for u in G.nodes():
        if 'label' not in G.node[u] or  G.node[u]['label'] not in ['synapse', 'root']:
            G.node[u]['label'] = 'steiner_midpoint'

    return G

def build_neuron(algorithm='snider', dim=3, alpha=0.5, **kwargs):
    #unmarked_points = grid_points3d(xmin=-2, xmax=2, ymin=-2, ymax=2, zmin=0, zmax=0)
    unmarked_points = random_points_uniform(num_points=400, xmin=-10, xmax=10, ymin=-10, ymax=10, zmin=0, zmax=0)
    G = init_graph(dim=dim)
    done = False
    while not done:
        if len(unmarked_points) == 0:
            break
        can_extend = update_graph(G, algorithm, unmarked_points, alpha, dim=dim, **kwargs)
        done = not can_extend
        
    retract_graph(G)

    if G.number_of_nodes() < 10:
        return None

    trees_dir = '%s/%s/trees' % (OUTDIR, algorithm)
    print trees_dir
    os.system('mkdir -p %s' % trees_dir)
    tree_dir = '%s/tree%d' % (trees_dir, len(os.listdir(trees_dir)) + 1)
    print tree_dir
    os.system('mkdir -p %s' % tree_dir)
    write_to_swc(G, outfile='%s/tree.swc' % tree_dir)
    with open('%s/synapses.txt' % tree_dir, 'w') as f:
        synapses = G.graph['synapses']
        for synapse in synapses:
            f.write('%d\n' % synapse)

    with open('%s/parameters.txt' % tree_dir, 'w') as f:
        f.write('algorithm %s\n' % algorithm)
        f.write('dimension %d\n' % dim)
        f.write('alpha %f\n' % alpha)
        for key, value in kwargs.iteritems():
            f.write('%s %s\n' % (key, value))

    G = read_tree(tree_dir)
    #print G.nodes()
    #print G.graph['synapses']
            
def build_neuron_video(algorithm='snider', dim=3, alpha=0.5, **kwargs):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    
    num_points = kwargs['npoints']
    xmin = kwargs['xmin']
    xmax = kwargs['xmax']
    ymin = kwargs['ymin']
    ymax = kwargs['ymax']
    zmin = kwargs['zmin']
    zmax = kwargs['zmax']
    
    init_points = random_points_uniform(num_points=400, xmin=xmin, xmax=xmax,\
                                                        ymin=ymin, ymax=ymax,\
                                                        zmin=zmin, zmax=zmax)
    for i in xrange(len(init_points)):
        init_points[i] = init_points[i][:dim]
    points = init_points[:]
    G = init_graph(dim=dim)
    graphs = [G.copy()]
    unmarked_points = [init_points]

    done = False
    max_steps = kwargs['max_steps']
    for i in xrange(max_steps):
        if len(points) == 0:
            break
        can_extend = update_graph(G, algorithm, points, alpha, dim=dim, **kwargs)
        graphs.append(G.copy())
        unmarked_points.append(points[:])
        if not can_extend:
            break
        
    retract_graph(G)
    graphs.append(G.copy())
    unmarked_points.append(points[:])

    def init():
        x, y, z = zip(*init_points)
        ax.scatter(x, y, zs=z, s=25)
        
    def init2d():
        x, y = zip(*init_points)
        plt.scatter(x, y, s=250)
        
    def redraw(frame):
        plt.clf()
        G = graphs[frame]
        for u in G.nodes():
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
        for u in G.nodes():
            coord = G.node[u]['coord']
            pos[u] = (coord[0], coord[1])
        viz_tree(G, save=False)

    ani = animation.FuncAnimation(fig, redraw2d, init_func=init, frames=len(graphs), \
                                  interval = 200)
    mywriter = animation.AVConvWriter()
    ani.save('neuron_builder/neuron_builder.mp4', writer=mywriter)
            
    plt.close()
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--alpha', default=0.5, type=float)
    parser.add_argument('--algorithm', default='snider')
    parser.add_argument('-d', '--dim', default=3, type=int)
    parser.add_argument('-rp', '--radius_puncta', default=1, type=float)
    parser.add_argument('-rr', '--radius_remove', default=1, type=float)
    parser.add_argument('-l', '--trial_length', default=1, type=float)
    parser.add_argument('-ra', '--radius_annihilation', default=1, type=float)
    parser.add_argument('-p', '--branching_prob', default=0.1, type=float)
    parser.add_argument('-v', '--video', action='store_true')
    parser.add_argument('-proot', '--prob_root', default=0.1, type=float)
    parser.add_argument('-m', '--max_steps', default=MAX_STEPS, type=int)
    parser.add_argument('-xmin', type=float, default=-10)
    parser.add_argument('-xmax', type=float, default=10)
    parser.add_argument('-ymin', type=float, default=-10)
    parser.add_argument('-ymax', type=float, default=10)
    parser.add_argument('-zmin', type=float, default=0)
    parser.add_argument('-zmax', type=float, default=0)
    parser.add_argument('-n', '--npoints', type=int, default=400)
    
    args = parser.parse_args()
    algorithm = args.algorithm
    dim = args.dim
    alpha = args.alpha
    video = args.video
    kwargs = vars(args)
    del kwargs['algorithm'], kwargs['alpha'], kwargs['dim']

    if video:
        build_neuron_video(algorithm, dim, alpha, **kwargs)
    else:
        build_neuron(algorithm, dim, alpha, **kwargs)
    #G = init_graph()
    #add_bifurcations(G, 1, dim=3, bifurcations=1, dist=0.775567698885)

if __name__ == '__main__':
    main()
