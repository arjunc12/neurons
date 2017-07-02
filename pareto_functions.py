import networkx as nx
from neuron_utils import point_dist, pareto_cost
from kruskal import kruskal, random_mst
from itertools import combinations
import numpy as np
from cost_functions import *
import pylab
from bisect import bisect_left, insort
from random import sample, choice
from collections import defaultdict

POP_SIZE = 50
GENERATIONS = 500

def satellite_tree(G):
    root = G.graph['root']
    
    satellite = G.copy()
    satellite.remove_edges_from(G.edges())

    for u in satellite.nodes_iter():
        if u != root:
            satellite.add_edge(u, root)
            p1, p2 = satellite.node[u]['coord'], satellite.node[root]['coord']
            satellite[u][root]['length'] = point_dist(p1, p2)
    
    return satellite

def min_spanning_tree(G):
    '''
    sorted_edges = sorted(G.edges(), key= lambda x : G[x[0]][x[1]]['length'])
    mst_edges = kruskal(G.nodes(), sorted_edges)
    mst = G.copy()
    mst.remove_edges_from(mst.edges())
    for u, v in mst_edges:
        mst.add_edge(u, v)
        mst[u][v]['length'] = G[u][v]['length']
    return mst
    '''
    return nx.minimum_spanning_tree(G, weight='length')

def replace_edge(G, T, u1, v1, u2, v2):
    T.remove_edge(u1, v1)
    T.add_edge(u2, v2)
    T[u2][v2]['length'] = G[u2][v2]['length']

def crossover_trees(G, p1, p2):
    T1 = p1.copy()
    T2 = p2.copy()

    u1, v1 = choice(T1.edges())
    u2, v2 = choice(T2.edges())

    replace_edge(G, T1, u1, v1, u2, v2)
    replace_edge(G, T2, u2, v2, u1, v1)

    return T1, T2

def pareto_genetic(G, alpha, axon=False, pop_size=POP_SIZE, generations=GENERATIONS):
    population = []
    for i in xrange(pop_size):
        mst = random_mst(G)
        mcost, scost = graph_costs(mst)
        cost = pareto_cost(mcost, scost, alpha)
        population.append((cost, mst))

    population = sorted(population)

    for generation in xrange(generations):
        parent1, parent2 = sample(population, 2)

        c1, p1 = parent1
        c2, p2 = parent2

        o1, o2 = crossover_trees(G, p1, p2)

        mcost1, scost1 = graph_costs(o1)
        cost1 = pareto_cost(mcost1, scost1, alpha)
        mcost2, scost2 = graph_costs(o2)
        cost2 = pareto_cost(mcost2, scost2, alpha)

        offspring1 = (cost1, o1)
        offspring2 = (cost2, o2)

        best_offspring = min([offspring1, offspring2])
       
        insort(population, best_offspring)
        
        worst_ind = population.pop()

    best_cost, best_tree = population[0]
    return best_tree

def pareto_kruskal(G, alpha, axon=False):
    root = G.graph['root']
    H = nx.Graph()
   
    H.add_node(root)
    H.graph['root'] = root
    H.node[root]['droot'] = 0
   
    graph_mcost = 0
    graph_scost = 0

    closest_neighbors = {}
    for u in G.nodes_iter():
        closest_neighbors[u] = G.node[u]['close_neighbors'][:]

    unpaired_neighbors = []
    candidate_nodes = defaultdict(int)

    while H.number_of_nodes() < G.number_of_nodes():
        best_edge = None
        best_cost = float("inf")

        #candidate_edges = set()
        candidate_edges = []
        for u in H.nodes():
            if axon and (u == H.graph['root']) and (H.degree(u) > 0):
                continue

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
                #candidate_edges.add(tuple(sorted((u, closest_neighbor))))
                candidate_edges.append((u, closest_neighbor))
                candidate_nodes[closest_neighbor].append(u)

        for u, v in candidate_edges:
            length = G[u][v]['length']
            
            #mcost = length + graph_mcost 
            mcost = length

            #scost = graph_scost
            scost = 0
            
            assert H.has_node(u)
            scost += length + H.node[u]['droot']
            '''
            if H.has_node(u):
                scost += length + H.node[u]['droot']
            elif H.has_node(v):
                scost += length + H.node[v]['droot'] 
            else:
                raise ValueError('something is wrong')
            '''

            
            cost = pareto_cost(mcost, scost, alpha)
            
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
        if best_edge == None:
            break
        u, v = best_edge
        H.add_edge(u, v)

        assert 'droot' in H.node[u]
        H.node[v]['droot'] = H.node[u]['droot'] + G[u][v]['length']
        
        '''
        if 'droot' in H.node[u]:
            H.node[v]['droot'] = H.node[u]['droot'] + G[u][v]['length']
        elif 'droot' in H.node[v]:
            H.node[u]['droot'] = H.node[v]['droot'] + G[u][v]['length']
        else:
            raise ValueError('something is really wrong')
        '''

        closest_neighbors[u].remove(v)
        closest_neighbors[v].append(u)
   
    pareto_mst = G.copy()
    pareto_mst.remove_edges_from(G.edges())
    for u, v in H.edges():
        pareto_mst.add_edge(u, v)
        pareto_mst[u][v]['length'] = G[u][v]['length']

    return pareto_mst

def initialize_khuller(mst):
    root = mst.graph['root']
    mst.node[root]['parent'] = None
    mst.node[root]['droot'] = 0
    queue = [root]
    while len(queue) > 0:
        curr = queue.pop(0)
        for child in mst.neighbors(curr):
            if child != mst.node[curr]['parent']:
                mst.node[child]['parent'] = curr
                queue.append(child)

def pareto_khuller(G, alpha, mst=None):
    if mst == None:
        mst = nx.minimum_spanning_tree(G, weight='length')
        
    initialize_khuller(mst)
    
    H = mst.copy()

    root = H.graph['root']
    queue = [root]
    #root_coord = H.node[root]['coord']
    while len(queue) > 0:
        curr = queue.pop()
        parent = H.node[curr]['parent']
        if parent != None:
            best_cost = float('inf')
            best_parent = None
            for candidate in nx.shortest_path(H, parent, root):
                wt = G[curr][candidate]['length']
                dist = wt + H.node[candidate]['droot']

                cost = pareto_cost(wt, dist, alpha)

                if cost < best_cost:
                    best_cost = cost
                    best_parent = candidate

            H.remove_edge(curr, parent)
            H.add_edge(curr, best_parent)
            H.node[curr]['parent'] = best_parent
            H[curr][best_parent]['length'] = G[curr][best_parent]['length']
            H.node[curr]['droot'] = H[curr][best_parent]['length'] + H.node[best_parent]['droot']

        for child in H.neighbors(curr):
            if H.node[child]['parent'] == curr:
                queue.append(child)
                H.node[child]['droot'] = H.node[curr]['droot'] + H[curr][child]['length']

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
        print alpha
        pareto_mst = pareto_khuller(G, T_span, alpha)
        print "satellite cost", satellite_cost(pareto_mst)
        print "mst cost", mst_cost(pareto_mst)

if __name__ == '__main__':
    main()
