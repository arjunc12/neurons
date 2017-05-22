import networkx as nx
from sys import argv
from itertools import combinations

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
        add_path(G, KT, T_sat, parent)
        relax(G, KT, parent, u)

def DFS(G, KT, T_span, T_sat, u, alpha):
    root = KT.graph['root']
    if KT.node[u]['droot'] > alpha * nx.shortest_path_length(T_sat, root, u, weight='length'):
        add_path(G, KT, T_sat, u)

    for v in T_span.neighbors(u):
        if v != T_span.node[u]['parent']:
            relax(G, KT, u, v)
            DFS(G, KT, T_span, T_sat, v, alpha)

def khuller(G, T_span, T_sat, root, alpha):
    assert alpha > 1
    initialize_parents(T_span)
    initialize_parents(T_sat)
    KT = T_span.copy()
    initialize(KT)
    DFS(G, KT, T_span, T_sat, G.graph['root'], alpha)
    KT.remove_edges_from(G.edges())
    for u in KT.nodes_iter():
        if u != KT.graph['root']:
            parent = KT.node[u]['parent']
            assert parent != None
            KT.add_edge(u, parent)

    return KT

def main():
    G = nx.Graph()
    root = (0, 0)
    points = [root, (0, 0.1), (0.001, 0.2), (0.2, 0), (1, 1), (2, 2.000001), (2, 2.1), (2.1, 2), (0.1001, 0.1)]
    for p1, p2 in combinations(points, 2):
        G.add_edge(p1, p2)
        G[p1][p2]['length'] = (((p1[0] - p2[0]) ** 2) + ((p1[1] - p2[1]) ** 2)) ** 0.5

    T_sat = nx.Graph()
    for u in G.nodes_iter():
        if u != root:
            T_sat.add_edge(u, root)
            T_sat[u][root]['length'] = (((u[0] - root[0]) ** 2) + ((u[1] - root[1]) ** 2)) ** 0.5

    T_span = nx.minimum_spanning_tree(G, weight='length')
    for H in [G, T_span, T_sat]:
        H.graph['root'] = root

    alpha = float(argv[1])

    KT = khuller(G, T_span, T_sat, root, alpha)
    print sorted(map(lambda x : tuple(sorted(x)), KT.edges()))


if __name__ == '__main__':
    main()
