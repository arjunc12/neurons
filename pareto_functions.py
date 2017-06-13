import networkx as nx

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
    H.node[root]['parent'] = None
    H.node[root]['droot'] = 0

def pareto_khuller(G, alpha, mst=None):
    if mst == None:
        mst = nx.minimum_spanning_tree(G, weight='length')
    H = mst.copy()
    

