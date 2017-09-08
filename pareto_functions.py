import networkx as nx
from neuron_utils import point_dist, pareto_cost, sort_neighbors
from kruskal import kruskal
from random_graphs import random_mst
from itertools import combinations
import numpy as np
from cost_functions import *
import pylab
from bisect import bisect_left, insort
from random import sample, choice, uniform, randint, random
from collections import defaultdict
from prufer import *
from enumerate_trees import find_all_spanning_trees
from random_graphs import random_point_graph

POP_SIZE = 400
GENERATIONS = 20000
MUTATION_PROB = 0.01

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

def crossover_trees(seq1, seq2):
    assert len(seq1) == len(seq2)
    
    s1 = seq1[:]
    s2 = seq2[:]
    
    i = randint(0, len(s1) - 1)
    j = randint(0, len(s2) - 1)

    s1[i] = seq2[j]
    s2[j] = seq1[i]

    return s1, s2

def mutate_tree(seq1, mutation_prob=MUTATION_PROB):
    if random() < MUTATION_PROB:
        idx = randint(0, len(seq1) - 1)
        new_digit = randint(0, 9)
        seq1[idx] = new_digit

def seq_to_tree(seq, G):
    T = from_prufer_sequence(seq)
    T.graph = G.graph
    for u in T.nodes():
        T.node[u] = G.node[u]

    for u, v in T.edges():
        T[u][v] = G[u][v]
        T[v][u] = G[v][u]
    
    return T

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

def pareto_genetic(G, axon=False, pop_size=POP_SIZE, generations=GENERATIONS,\
                   mutation_prob=MUTATION_PROB):
    root = G.graph['root']
    
    label_map = {}
    label_map_inv = {}
    for i, u in enumerate(G.nodes()):
        label_map[u] = i
        label_map_inv[i] = u

    G = nx.relabel_nodes(G, label_map)
    G.graph['root'] = label_map[root]

    population = []
    for i in xrange(pop_size):
        mst = random_mst(G)
        mcost, scost = graph_costs(mst)
        #cost = pareto_cost(mcost, scost, alpha)
        population.append((mst, mcost, scost))

    #population = sorted(population)

    for generation in xrange(generations):
        for i, (ind, mcost, scost) in enumerate(population):
            if random() < MUTATION_PROB:
                seq = to_prufer_sequence(ind)
                mutate_tree(seq, mutation_prob=1)
                ind = seq_to_tree(seq, G)
                population[i] = (ind, mcost, scost)

        parent1, parent2 = sample(population, 2)

        p1, m1, s1 = parent1
        p2, m2, s2 = parent2

        p1 = to_prufer_sequence(p1)
        p2 = to_prufer_sequence(p2)
        #mutate_tree(p1, MUTATION_PROB)
        #mutate_tree(p2, MUTATION_PROB)

        o1, o2 = crossover_trees(p1, p2)
        #mutate_tree(o1, MUTATION_PROB)
        #mutate_tree(o2, MUTATION_PROB)

        o1 = seq_to_tree(o1, G)
        o2 = seq_to_tree(o2, G)

        mcost1, scost1 = graph_costs(o1)
        #cost1 = pareto_cost(mcost1, scost1, alpha)
        mcost2, scost2 = graph_costs(o2)
        #cost2 = pareto_cost(mcost2, scost2, alpha)
        offspring = [(o1, mcost1, scost1), (o2, mcost2, scost2)]
        origin = (0, 0)
        for ospring, mcosto, scosto in offspring:
            worst_index = None
            worst_dist = 0
            for i, (ind, mcosti, scosti) in enumerate(population):
                if partially_dominates((mcosto, scosto), (mcosti, scosti)):
                    worst_index = i
                    break
            if worst_index != None:
                population[worst_index] = (ospring, mcosto, scosto)



        #offspring1 = (cost1, o1)
        #offspring2 = (cost2, o2)

        #best_offspring = min([offspring1, offspring2])
       
        #insort(population, best_offspring)
        
        #worst_ind = population.pop()

    best_trees = []
    for ind, mcost, scost in population:
        T = nx.relabel_nodes(ind, label_map_inv)
        T.graph['root'] = root
        best_trees.append((mcost, scost, T))
    return best_trees
    
    #best_cost, best_tree = population[0]
    #best_tree = nx.relabel_nodes(best_tree, label_map_inv)
    #best_tree.graph['root'] = root
    #return best_tree

def pareto_prim(G, alpha, axon=False):
    root = G.graph['root']
    if 'sorted' not in G.graph:
        sort_neighbors(G)

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
    candidate_nodes = defaultdict(list)

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

def pareto_prim_sandbox(G, alpha, axon=False):
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

    unpaired_nodes = [root]
    candidate_nodes = defaultdict(list)

    while H.number_of_nodes() < G.number_of_nodes():
        best_edge = None
        best_cost = float("inf")

        candidate_edges = []
        for u in unpaired_nodes:
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
                #candidate_edges.append((u, closest_neighbor))
                candidate_nodes[closest_neighbor].append(u)

                length = G[u][closest_neighbor]['length']
                mcost = length
                scost = 0
                assert H.has_node(u)
                scost = length + H.node[u]['droot']
                cost = pareto_cost(mcost, scost, alpha)
                insort(candidate_edges, (cost, (u, closest_neighbor)))

        '''
        for v in candidate_nodes:
            for u in candidate_nodes[v]:
                length = G[u][v]['length']
                
                mcost = length

                scost = 0
                
                assert H.has_node(u)
                scost += length + H.node[u]['droot']
                
                cost = pareto_cost(mcost, scost, alpha)
                
                if cost < best_cost:
                    best_edge = (u, v)
                    best_cost = cost
        '''
        #print candidate_edges
        best_cost, best_edge = candidate_edges.pop(0)

        if best_edge == None:
            break
        u, v = best_edge
        H.add_edge(u, v)

        unpaired_nodes = [v] + candidate_nodes[v]
        del candidate_nodes[v]

        assert 'droot' in H.node[u]
        H.node[v]['droot'] = H.node[u]['droot'] + G[u][v]['length']
        
        closest_neighbors[u].remove(v)
        closest_neighbors[v].append(u)
   
    pareto_mst = G.copy()
    pareto_mst.remove_edges_from(G.edges())
    for u, v in H.edges():
        pareto_mst.add_edge(u, v)
        pareto_mst[u][v]['length'] = G[u][v]['length']

    return pareto_mst

def alpha_to_beta(alpha, opt_mcost, opt_scost):
    assert alpha > 0
    assert alpha < 1
    num = 2 * alpha * opt_mcost
    denom = (1 - alpha) * opt_scost
    frac = num / denom
    frac **= 0.5
    return 1 + frac

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

def pareto_brute_force(G, alpha, trees=None):
    if trees == None:
        trees = find_all_spanning_trees(G)

    best_tree = None
    best_cost = float("inf")
    for tree in trees:
        mcost, scost = graph_costs(tree)
        cost = pareto_cost(mcost, scost, alpha)
        if cost < best_cost:
            best_cost = cost
            best_tree = tree

    return best_tree

def main():
    G = random_point_graph(15)
    
    sat_tree = satellite_tree(G)
    span_tree = nx.minimum_spanning_tree(G, weight='length')
           
    opt_scost = float(satellite_cost(sat_tree))
    normalize_scost = lambda cost : normalize_cost(cost, opt_scost)
    opt_mcost = float(mst_cost(span_tree))
    normalize_mcost = lambda cost : normalize_cost(cost, opt_mcost)

    mcosts1  = []
    scosts1 = []

    sort_neighbors(G)

    delta = 0.01
    alphas = pylab.arange(delta, 1, delta)
    for alpha in alphas:
        print alpha
        pareto_mst = pareto_prim(G, alpha)
        mcost1, scost1 = graph_costs(pareto_mst)
        mcost1 = normalize_mcost(mcost1)
        scost1 = normalize_scost(scost1)
        mcosts1.append(mcost1)
        scosts1.append(scost1)

    best_trees = pareto_genetic(G)
    mcosts2 = []
    scosts2 = []
    for mcost2, scost2, tree in best_trees:
        mcost2 = normalize_mcost(mcost2)
        scost2 = normalize_scost(scost2)
        mcosts2.append(mcost2)
        scosts2.append(scost2)

    pylab.figure()
    pylab.plot(mcosts1, scosts1, c='b')
    pylab.scatter(mcosts1, scosts1, c='b')
    pylab.scatter(mcosts2, scosts2, c='r')
    pylab.savefig('genetic.pdf')
    pylab.close()

if __name__ == '__main__':
    main()
