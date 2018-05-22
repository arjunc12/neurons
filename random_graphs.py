import networkx as nx
from random import uniform, choice
from itertools import combinations
from neuron_utils import point_dist
from collections import defaultdict
from random import choice, shuffle
from numpy.random import permutation
from graph_utils import is_tree
from dist_functions import node_dist

def random_points(num_points=10, xmin=-10, xmax=10, ymin=-10, ymax=10,\
                                 zmin=-10, zmax=10):
    points = []
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
    G.graph['root'] = choice(G.nodes())
    
    for u, v in combinations(G.nodes(), 2):
        G.add_edge(u, v)
        G[u][v]['length'] = point_dist(points[u], points[v])

    return G

def random_mst(G, euclidean=False):
    H = G.copy()
    H.remove_edges_from(G.edges())
    root = H.graph['root']

    in_tree = defaultdict(bool)
    in_tree[root] = True
    successor = {}
    successor[root] = None

    vertices = H.nodes()
    shuffle(vertices)
    for u in vertices:
        curr = u
        while not in_tree[curr]:
            candidates = None
            if euclidean:
                candidates = G.nodes()
                candidates.remove(u)
            else:
                candidates = G.neighbors(curr)
            successor[curr] = choice(candidates)
            curr = successor[curr]

        curr = u
        while not in_tree[curr]:
            in_tree[curr] = True
            next = successor[curr]
            H.add_edge(curr, next)
            if euclidean:
                H[curr][next]['length'] = node_dist(G, curr, next)
            else:
                H[curr][next]['length'] = G[curr][next]['length']
            curr = next
    
    return H

def barabasi_tree(G):
    ba = nx.barabasi_albert_graph(n=G.number_of_nodes(), m=1)
    node_map = {}
    for i, u in enumerate(permutation(G.nodes())):
        node_map[i] = u

    H = G.copy()
    H.remove_edges_from(G.edges())

    for a, b in ba.edges_iter():
        u, v = node_map[a], node_map[b]
        H.add_edge(u, v)
        H[u][v]['length'] = node_dist(H, u, v)

    return H

def main():
    G = random_point_graph()
    G.remove_edges_from(G.edges())
    T = barabasi_tree(G)
    T2 = random_mst(G, euclidean=True)
    assert is_tree(T)
    assert is_tree(T2)

if __name__ == '__main__':
    main()
