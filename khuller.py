import networkx as nx
from sys import argv
from itertools import combinations
import numpy as np
from cost_functions import *
from neuron_utils import *

TEST_NEURON = 'datasets/bipolar/rat/neocortex/markram/RP140319-CHC-3-IDA.CNG.swc'

def initialize(KT):
    root = KT.graph['root']
    for u in KT.nodes_iter():
        if u == KT.graph['root']:
            KT.node[u]['droot'] = 0
        else:
            KT.node[u]['droot'] = float('inf')
        if u == root:
            KT.node[u]['parent'] = None
        else:
            KT.node[u]['parent'] = root

def relax(G, KT, u, v):
    d1 = KT.node[v]['droot']
    d2 = KT.node[u]['droot'] + G[u][v]['length']
    if d1 > d2:
        KT.node[v]['droot'] = d2
        KT.node[v]['parent'] = u

def initialize_parents(T):
    root = T.graph['root']
    queue = [root]
    T.node[root]['parent'] = None
    visited = set()
    while len(queue) > 0:
        curr = queue.pop(0)
        assert curr not in visited
        visited.add(curr)
        for child in T.neighbors(curr):
            if T.node[curr]['parent'] == child:
                assert child in visited
            else:
                assert child not in visited
                T.node[child]['parent'] = curr
                queue.append(child)
    assert len(visited) == T.number_of_nodes()

def add_path(G, KT, T_sat, u):
    root = G.graph['root']
    if KT.node[u]['droot'] > nx.shortest_path_length(T_sat, root, u, weight='length'):
        parent = T_sat.node[u]['parent']
        assert T_sat.has_node(parent)
        add_path(G, KT, T_sat, parent)
        relax(G, KT, parent, u)

def DFS(G, KT, T_span, T_sat, u, beta):
    root = KT.graph['root']
    if KT.node[u]['droot'] > beta * nx.shortest_path_length(T_sat, root, u, weight='length'):
        add_path(G, KT, T_sat, u)

    for v in T_span.neighbors(u):
        if v != T_span.node[u]['parent']:
            relax(G, KT, u, v)
            DFS(G, KT, T_span, T_sat, v, beta)


def khuller(G, T_span, T_sat, beta):
    assert beta > 1
    initialize_parents(T_span)
    initialize_parents(T_sat)
    KT = T_span.copy()
    initialize(KT)
    DFS(G, KT, T_span, T_sat, G.graph['root'], beta)
    KT.remove_edges_from(G.edges())
    for u in KT.nodes_iter():
        if u != KT.graph['root']:
            parent = KT.node[u]['parent']
            assert parent != None
            KT.add_edge(u, parent)
            KT[u][parent]['length'] = G[u][parent]['length']
    
    assert nx.is_connected(KT)
    assert KT.number_of_edges() == KT.number_of_nodes() - 1
    return KT

def khuller_test_graph():
    G = nx.Graph()
    root = (0, 0)
    G.add_node(root)
    G.graph['root'] = root
    G.node[root]['coord'] = (0, 0)

    for i in xrange(1, 11):
        coord1 = (1, 2 + i)
        coord2 = (-1, 2 + i)
        G.add_node(coord1)
        G.node[coord1]['coord'] = coord1
        G.add_node(coord2)
        G.node[coord2]['coord'] = (-1, 2 + i)

    return complete_graph(G)

def main():
    G = khuller_test_graph()
    T_sat = satellite_tree(G)
    opt_scost = satellite_cost(T_sat)
    T_span = min_spanning_tree(G)
    opt_mcost = mst_cost(T_span)

    alpha = float(argv[1])
    beta = alpha_to_beta(alpha, opt_mcost, opt_scost)
    print beta

    KT = khuller(G, T_span, T_sat, beta)
    for u, v in sorted(KT.edges_iter()):
        print KT.node[u]['coord'], KT.node[v]['coord']

if __name__ == '__main__':
    main()
