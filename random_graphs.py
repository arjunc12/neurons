import networkx as nx
from random import uniform
from itertools import combinations
from neuron_utils import point_dist
from collections import defaultdict
from random import choice

def random_points(num_points=10, xmin=-10, xmax=10, ymin=-10, ymax=10,\
                                 zmin=-10, zmax=10):
    points = [(0, 0, 0)]
    for i in xrange(num_points):
        x = uniform(xmin, xmax)
        y = uniform(ymin, ymax)
        z = uniform(zmin, zmax)
        points.append((x, y, z))

    return points

def random_point_graph(num_points=10, xmin=-10, xmax=10,\
                                      ymin=-10, ymax=10,\
                                      zmin=-10, zmax=10):
    points = random_points(num_points, xmin, xmax, ymin, ymax, zmin, zmax)
    G = nx.Graph()
    for i, point in enumerate(points):
        G.add_node(i)
        G.node[i]['coord'] = point
    G.graph['root'] = 0
    
    for u, v in combinations(G.nodes(), 2):
        G.add_edge(u, v)
        G[u][v]['length'] = point_dist(points[u], points[v])

    return G

def random_mst(G):
    H = G.copy()
    H.remove_edges_from(G.edges())
    root = H.graph['root']

    in_tree = defaultdict(bool)
    in_tree[root] = True
    successor = {}
    successor[root] = None

    for u in H.nodes():
        curr = u
        while not in_tree[curr]:
            successor[curr] = choice(G.neighbors(curr))
            curr = successor[curr]

        curr = u
        while not in_tree[curr]:
            in_tree[curr] = True
            next = successor[curr]
            H.add_edge(curr, next)
            H[curr][next]['length'] = G[curr][next]['length']
            curr = next
    
    return H
