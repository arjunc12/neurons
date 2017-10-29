import networkx as nx
from neuron_utils import point_dist

def partially_dominates(costs1, costs2):
    assert len(costs1) == len(costs2)
    strict = False
    for i in xrange(len(costs1)):
        cost1, cost2 = costs1[i], costs2[i]
        if cost1 > cost2:
            return False
        elif cost1 < cost2:
            strict = True

    return strict

def prop_dominated(mcosts1, scosts1, mcosts2, scosts2):
    assert len(mcosts1) == len(scosts1)
    assert len(mcosts2) == len(scosts2)
    comparisons = 0
    dominated = 0
    for i in xrange(len(mcosts2)):
        comparisons += 1
        costs2 = (mcosts2[i], scosts2[i])
        for j in xrange(len(mcosts1)):
            costs1 = (mcosts1[j], scosts1[j])
            if partially_dominates(costs1, costs2):
                dominated += 1
                break
    return comparisons, dominated

def graph_normalize_functions(G):
    opt_mst_cost = best_mst_cost(G)
    opt_sat_cost = best_satellite_cost(G)

    mst_func = make_normalize_function(opt_mst_cost)
    sat_func = make_normalize_function(opt_sat_cost)

    return mst_func, sat_func

def make_normalize_function(opt_cost):
    return lambda cost : normalize_cost(cost, opt_cost)

def normalize_cost(cost, opt_cost):
    return 1 - (opt_cost / cost)

def graph_costs(G):
    scost = 0
    mcost = 0
    
    droot = {}
    root = G.graph['root']
    droot[root] = 0
    
    parent = {}
    parent[root] = None

    queue = [root]
    visited = set()

    mcosts = []
    scosts = []
    while len(queue) > 0:
        curr = queue.pop(0)
        
        if curr in visited:
            return float("inf"), float("inf")
        
        visited.add(curr)
        for child in G.neighbors(curr):
            if child != parent[curr]:
                length = G[curr][child]['length']
                #mcost += length
                mcosts.append(length)
                child_droot = length + droot[curr]
                droot[child] = child_droot
                #scost += child_droot
                scosts.append(child_droot)
                parent[child] = curr
                queue.append(child)

    if len(visited) < G.number_of_nodes():
        scost = float("inf") # graph is not connected

    assert len(visited) == G.number_of_nodes()

    mcost = sum(sorted(mcosts))
    scost = sum(sorted(scosts))

    return mcost, scost

def best_mst_cost(G):
    mst = nx.minimum_spanning_tree(G, weight='length')
    best_cost = 0
    for u, v in mst.edges_iter():
        best_cost += G[u][v]['length']
    return best_cost

def mst_cost(G):
    total_length = 0
    for u, v in G.edges_iter():
        total_length += G[u][v]['length']
    return total_length

def best_satellite_cost(G):
    best_cost = 0
    root = G.graph['root']
    root_coord = G.node[root]['coord']
    costs = []
    for u in G.nodes():
        if u == root:
            continue
        cost = point_dist(root_coord, G.node[u]['coord'])
        #best_cost += cost
        costs.append(cost)
    #return best_cost
    return sum(sorted(costs))
    #return KahanSum(costs)

def satellite_cost(G, relevant_nodes=None):
    root = G.graph['root']
    total_cost = 0
    costs = []
    if relevant_nodes == None:
        relevant_nodes = G.nodes()
    for u in relevant_nodes:
        if not G.has_node(u):
            continue
        if u == root:
            continue
        if 'droot' in G.node[u]:
            costs.append(G.node[u]['droot'])
        elif nx.has_path(G, root, u):
            droot = nx.shortest_path_length(G, root, u, weight='length')
            #total_cost += droot
            costs.append(droot)
        else:
            #total_cost += float("inf")
            costs.append(float("inf"))
    
    #return total_cost
    return sum(sorted(costs))
    #return KahanSum(costs)

def centroid_mst_costs(G):
    root = G.graph['root']
    centroidp = centroid(G)
    root_cdist = point_dist(root, centroidp)

    mcost = root_cdist
    scost = 0
    
    for point in points:
        if point != root:
            cdist = point_dist(point, centroidp)
            mcost += cdist
            scost += cdist + root_cdist

    return mcost, scost

def main():
    mcosts1 = [1, 1, 0, 100]
    scosts1 = [3, 4, 0, 100]

    mcosts2 = [0, 2]
    scosts2 = [2, 10]

    print prop_dominated(mcosts1, scosts1, mcosts2, scosts2)
    print prop_dominated(mcosts2, scosts2, mcosts1, scosts1)

if __name__ == '__main__':
    main()
