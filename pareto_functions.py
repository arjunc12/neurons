import networkx as nx
from neuron_utils import point_dist, pareto_cost

def pareto_kruskal(G, alpha):
    root = G.graph['root']
    H = nx.Graph()
   
    H.add_node(root)
    H.node[root]['droot'] = 0
   
    graph_mcost = 0
    graph_scost = 0

    closest_neighbors = {}
    for u in G.nodes_iter():
        closest_neighbors[u] = G.node[u]['close_neighbors'][:]

    while H.number_of_nodes() < G.number_of_nodes():
        best_edge = None
        best_cost = float("inf")

        candidate_edges = set()
        for u in H.nodes():
            assert 'droot' in H.node[u]
             
            invalid_neighbors = []
            closest_neighbor = None
            for i in xrange(len(closest_neighbors[u])):
                v = closest_neighbors[u][i]
                if H.has_node(v):
                    invalid_neighbors.append(v)
                else:
                    closest_neighbor = v
                    break

            for n in invalid_neighbors:
                closest_neighbors[u].remove(n)

            if closest_neighbor != None:
                candidate_edges.add(tuple(sorted((u, closest_neighbor))))

        for u, v in candidate_edges:

            scost = graph_scost           
            if H.has_node(u):
                scost += H.node[u]['droot'] + G[u][v]['length']
            elif H.has_node(v):
                scost += H.node[v]['droot'] + G[u][v]['length']
            else:
                raise ValueError('something is wrong')

            mcost = graph_mcost + G[u][v]['length']
            
            cost = pareto_cost(mcost, scost, alpha)
            
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
        if best_edge == None:
            break
        u, v = best_edge
        H.add_edge(u, v)

        if 'droot' in H.node[u]:
            H.node[v]['droot'] = H.node[u]['droot'] + G[u][v]['length']
        elif 'droot' in H.node[v]:
            H.node[u]['droot'] = H.node[v]['droot'] + G[u][v]['length']
        else:
            raise ValueError('something is really wrong')

        closest_neighbors[u].remove(v)
        closest_neighbors[v].append(u)
   
    pareto_mst = G.copy()
    pareto_mst.remove_edges_from(G.edges())
    for u, v in H.edges():
        pareto_mst.add_edge(u, v)
        pareto_mst[u][v]['length'] = G[u][v]['length']

    return pareto_mst

def initialize_khuller(mst):
    root = H.graph['root']
    mst.node[root]['parent'] = None
    mst.node[root]['droot'] = 0
    queue = [root]
    while len(queue) > 0:
        curr = queue.pop(0)
        for child in mst.neighbors(curr):
            if child != mst.node[curr]['parent']:
                mst.node[child]['parent'] = curr
                queue.append(node)

def pareto_khuller(G, alpha, mst=None):
    if mst == None:
        mst = nx.minimum_spanning_tree(G, weight='length')
        initialize_khuller(mst)
    
    H = mst.copy()

    root = H.graph['root']
    queue = [root]
    root_coord = H.node[root]['coord']
    while len(queue) > 0:
        curr = queue.pop()
        parent = H.node[curr]['parent']
        if parent != None:
            best_cost = float("inf")
            best_parent = None
            for candidate in nx.shortest_path(H, parent, root):
                w1 = H[curr][candidate]['next']
                w2 = G.node[curr][root]['length']

                d1 = H.node[curr]['droot']
                d2 = w2

                cost1 = pareto_cost(w1, d1, alpha)
                cost2 = pareto_cost(w1, d2, alpha)

                if cost2 < cost1:
                    H.node[curr]['parent'] = root
                    H.node[curr]['droot']
                    H[curr][root]['length'] = w2
                    H.remove_edge(curr, parent)

        for child in H.neighbors(curr):
            if H.node[child]['parent'] == curr:
                queue.append(child)
                H.node[child]['droot'] = H.node[curr]['droot'] + H[curr][child]['weight']

    return H

def main():
    pass
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
    
    for alpha in np.arange(0.01, 0.99, 0.01):
        print beta
        pareto_mst = pareto_khuller(G, alpha, T_span)
        print "satellite cost", satellite_cost(KT)
        print "mst cost", mst_cost(KT)

if __name__ == '__main__':
    pass
