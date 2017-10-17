import networkx as nx
from neuron_utils import *
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
from steiner_midpoint import *
from sys import argv
from scipy.spatial.distance import euclidean
from graph_utils import *

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
        new_digit = randint(0, len(seq1))
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

    best_trees = []
    for ind, mcost, scost in population:
        T = nx.relabel_nodes(ind, label_map_inv)
        T.graph['root'] = root
        best_trees.append((mcost, scost, T))
    return best_trees    

def pareto_prim(G, alpha, axon=False):
    root = G.graph['root']
    if 'sorted' not in G.graph:
        sort_neighbors(G)

    H = nx.Graph()
   
    H.add_node(root)
    H.graph['root'] = root
    H.node[root]['droot'] = 0
    H.node[root]['coord'] = G.node[root]['coord']
   
    graph_mcost = 0
    graph_scost = 0

    closest_neighbors = {}
    for u in G.nodes_iter():
        closest_neighbors[u] = G.node[u]['close_neighbors'][:]

    unpaired_neighbors = []
    candidate_nodes = defaultdict(list)

    while H.number_of_nodes() < G.number_of_nodes():
        best_edge = None
        best_mcost = None
        best_scost = None
        best_cost = float("inf")

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
                candidate_edges.append((u, closest_neighbor))
                candidate_nodes[closest_neighbor].append(u)

        for u, v in candidate_edges:
            length = G[u][v]['length']
            
            mcost = length

            scost = 0
            
            assert H.has_node(u)
            scost += length + H.node[u]['droot']
            
            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
            
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
                best_mcost = mcost
                best_scost = scost

        if best_edge == None:
            break
        u, v = best_edge

        H.add_edge(u, v)
        H[u][v]['length'] = G[u][v]['length']

        assert 'droot' in H.node[u]
        H.node[v]['droot'] = H.node[u]['droot'] + G[u][v]['length']
        H.node[v]['coord'] = G.node[v]['coord']

        check_dists(H)
 
        closest_neighbors[u].remove(v)
        closest_neighbors[v].remove(u)
         
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
        best_mcost = None
        best_scost = None
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
        closest_neighbors[v].remove(u)
        
 
    pareto_mst = G.copy()
    pareto_mst.remove_edges_from(G.edges())
    for u, v in H.edges():
        pareto_mst.add_edge(u, v)
        pareto_mst[u][v]['length'] = G[u][v]['length']

    return pareto_mst

def add_midpoint(G, u, parent, midpoint, midpoint_node):
    assert G.has_edge(u, parent)
    G.remove_edge(u, parent)
    G.add_node(midpoint_node)
    G.node[midpoint_node]['coord'] = midpoint
    G.node[midpoint_node]['label'] = 'steiner_point'

    G.add_edge(u, midpoint_node)
    G[u][midpoint_node]['length'] = node_dist(G, u, midpoint_node)
    G.node[u]['parent'] = midpoint_node

    G.add_edge(midpoint_node, parent)
    G[midpoint_node][parent]['length'] = node_dist(G, midpoint_node, parent)
    G.node[midpoint_node]['parent'] = parent

    G.node[midpoint_node]['droot'] = G[midpoint_node][parent]['length'] + G.node[parent]['droot']

def pareto_steiner(G, alpha, axon=False):
    root = G.graph['root']
    if 'sorted' not in G.graph:
        sort_neighbors(G)

    H = nx.Graph()
   
    H.add_node(root)
    H.graph['root'] = root
    H.node[root]['droot'] = 0
    H.node[root]['parent'] = None
    root_coord = G.node[root]['coord']
    H.node[root]['coord'] = root_coord
    added_nodes = 1

    in_nodes = set()
    out_nodes = set(G.nodes())
    in_nodes.add(root)
    out_nodes.remove(root)
   
    graph_mcost = 0
    graph_scost = 0

    closest_neighbors = {}
    for u in G.nodes_iter():
        closest_neighbors[u] = G.node[u]['close_neighbors'][:]

    unpaired_neighbors = []
    candidate_nodes = defaultdict(list)

    node_index = max(G.nodes()) + 1

    dist_error = 0

    steps = 0
    midpoints = {}
    choices = {}
    while added_nodes < G.number_of_nodes():
        best_edge = None
        best_mcost = None
        best_scost = None
        best_cost = float("inf")

        best_choice = None
        best_midpoint = None

        candidate_edges = []

        for u in in_nodes:
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
                candidate_edges.append((u, closest_neighbor))
                #candidate_edges.append((root, closest_neighbor))
                candidate_nodes[closest_neighbor].append(u)
 
        for u, v in candidate_edges:
            assert H.has_node(u)
            assert not H.has_node(v)

            p1 = H.node[u]['coord']
            p3 = G.node[v]['coord']
            midpoint = None
            choice = None
            w = H.node[u]['parent']
            if w != None:
                assert H.has_node(w)
                p2 = H.node[w]['coord']

            if (u, v) in midpoints:
                assert (u, v) in choices
                midpoint = midpoints[(u, v)]
                choice = choices[(u, v)]
            else:
                if w == None:
                    midpoint = p1
                    choice = 1
                else:
                    assert p2 != None
                    droot = H.node[w]['droot']
                    midpoint, choice = best_midpoint_approx(p1, p2, p3, alpha)
                assert midpoint != None
                assert choice != None
                midpoints[(u, v)] = midpoint
                choices[(u, v)] = choice
                
            length = point_dist(p3, midpoint)
            mcost = length
            scost = length
            if choice == 1:
                scost += H.node[u]['droot']
            else:
                assert w != None
                assert p2 != None
                scost += H.node[w]['droot']
                if choice == 3:
                    assert p2 != None
                    scost += point_dist(midpoint, p2)
           

            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
            
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
                best_midpoint = midpoint
                best_choice = choice
                best_mcost = mcost
                best_scost = scost
        

        if best_edge == None:
            break
        u, v = best_edge
        for a, b in midpoints.keys():
            assert (a, b) in choices
            if b == v:
                del midpoints[(a, b)]
                del choices[(a, b)]
        w = H.node[u]['parent']

        out_node = v
        in_node = None
        steiner = False
        if best_choice == 1:
            in_node = u
        elif best_choice == 2:
            in_node = H.node[u]['parent']
        else:
            assert best_choice == 3
            assert w != None
           

            p1 = H.node[u]['coord']
            p2 = H.node[w]['coord']

            steiner = True
            midpoint_node = node_index
            in_node = midpoint_node
            node_index += 1

            init_dist = H.node[u]['droot']
            init_length = H[u][w]['length']

            H.add_node(midpoint_node)
            H.node[midpoint_node]['coord'] = best_midpoint
            H.node[midpoint_node]['label'] = 'steiner_point'

            assert H.has_edge(u, w)
            H.remove_edge(u, w) 
            H.add_edge(midpoint_node, w)
            H.add_edge(u, midpoint_node)
            
            H[midpoint_node][w]['length'] = point_dist(best_midpoint, p2)
            H[u][midpoint_node]['length'] = point_dist(p1, best_midpoint)

            H.node[u]['parent'] = midpoint_node
            H.node[midpoint_node]['parent'] = w

            midpoint_dist = H[midpoint_node][w]['length'] + H.node[w]['droot']
            H.node[midpoint_node]['droot'] = midpoint_dist

            H.node[u]['droot'] = H[u][midpoint_node]['length'] + H.node[midpoint_node]['droot']

            midpoint_neighbors = []
            for u in G.nodes():
                if not H.has_node(u):
                    p1 = H.node[midpoint_node]['coord']
                    p2 = G.node[u]['coord']
                    dist = point_dist(p1, p2)
                    midpoint_neighbors.append((dist, u))

            midpoint_neighbors = sorted(midpoint_neighbors)
            closest_neighbors[midpoint_node] = []
            for dist, neighbor in midpoint_neighbors:
                closest_neighbors[midpoint_node].append(neighbor)

        H.add_node(out_node)
        H.node[out_node]['coord'] = G.node[out_node]['coord']
        H.add_edge(in_node, out_node)
        H.node[out_node]['parent'] = in_node
        in_coord = H.node[in_node]['coord']
        out_coord = H.node[out_node]['coord']
        H[in_node][out_node]['length'] = point_dist(in_coord, out_coord)

        assert 'droot' in H.node[in_node]
        H.node[out_node]['droot'] = H.node[in_node]['droot'] + H[in_node][out_node]['length']
       
        if True:
            if out_node in closest_neighbors[in_node]:
                closest_neighbors[in_node].remove(out_node)
            if in_node in closest_neighbors[out_node]:
                closest_neighbors[out_node].remove(in_node)

        added_nodes += 1
        in_nodes.add(out_node)
        out_nodes.remove(out_node)
        assert len(in_nodes) == added_nodes
        assert len(in_nodes) + len(out_nodes) == G.number_of_nodes()
        for n in in_nodes:
            assert G.has_node(n)
            assert H.has_node(n)


    assert is_tree(H)

    pareto_mst = G.copy()
    pareto_mst.remove_edges_from(G.edges())
    for u, v in H.edges():
        pareto_mst.add_edge(u, v)
        pareto_mst[u][v]['length'] = H[u][v]['length']

    return H

def pareto_steiner_sandbox(G, alpha, axon=False):
    root = G.graph['root']
    if 'sorted' not in G.graph:
        sort_neighbors(G)

    H = nx.Graph()
   
    H.add_node(root)
    H.graph['root'] = root
    H.node[root]['droot'] = 0
    H.node[root]['coord'] = G.node[root]['coord']
    H.node[root]['parent'] = None
    added_nodes = 1
 
    graph_mcost = 0
    graph_scost = 0

    closest_neighbors = {}
    for u in G.nodes_iter():
        closest_neighbors[u] = G.node[u]['close_neighbors'][:]

    unpaired_neighbors = []
    candidate_nodes = defaultdict(list)

    node_index = max(G.nodes()) + 1

    in_nodes = set()
    in_nodes.add(root)
    out_nodes = set(G.nodes())
    out_nodes.remove(root)

    while added_nodes < G.number_of_nodes():
        best_edge = None
        best_mcost = None
        best_scost = None
        best_cost = float("inf")

        #candidate_edges = set()
        candidate_edges = []
        for u in in_nodes:
            if axon and (u == H.graph['root']) and (H.degree(u) > 0):
                continue

            if not G.has_node(u):
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
                candidate_edges.append((u, closest_neighbor))
                candidate_nodes[closest_neighbor].append(u)

        for u, v in candidate_edges:
            assert H.has_node(u)
            assert not H.has_node(v)
            length = G[u][v]['length']
            
            mcost = length
            assert H.has_node(u)
            scost = length + H.node[u]['droot']
            
            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)

            if added_nodes == 440:
                if u == 194:
                    print mcost, scost, cost, best_cost
            
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
                best_mcost = mcost
                best_scost = scost

        if best_edge == None:
            break
        u, v = best_edge
        out_node = v
        child = u
        parent = H.node[child]['parent']
        in_node = u
        best_midpoint = None
        best_midpoint_node = None
        best_midpoint_edge = (child, parent)
        best_choice = 1
        p3 = G.node[v]['coord']

        #while parent != None:
        if parent != None:
            p1 = H.node[child]['coord']
            p2 = H.node[parent]['coord']
            
            midpoint, choice = best_midpoint_approx(p1, p2, p3, alpha)
            mcost = point_dist(p3, midpoint)
            scost = mcost + point_dist(midpoint, p2) + H.node[parent]['droot']

            if choice == 3 and scost > best_scost:
                print "weird midpoint scost"
                print best_scost, scost

            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)

            if cost < best_cost:
                best_mcost = mcost
                best_scost = scost
                best_cost = cost
                best_midpoint = midpoint
                best_choice = choice
                best_midpoint_edge = (child, parent)

            #child = parent
            #parent = H.node[parent]['parent']

        #assert child == root
        #assert parent == None

        child, parent = best_midpoint_edge
        if best_choice == 1:
            in_node = child
        elif best_choice == 2:
            in_node = parent
        else:
            assert best_choice == 3
            assert best_midpoint != None

            midpoint_node = node_index
            in_node = midpoint_node
            node_index += 1
            add_midpoint(H, child, parent, best_midpoint, midpoint_node)
        
        assert H.has_node(in_node)
        assert not H.has_node(out_node)
        assert G.has_node(out_node)

        H.add_edge(in_node, out_node)
        H.node[out_node]['coord'] = G.node[out_node]['coord']
        H[in_node][out_node]['length'] = node_dist(H, in_node, out_node)

        assert 'droot' in H.node[u]
        H.node[out_node]['droot'] = H.node[in_node]['droot'] + node_dist(H, in_node, out_node)
        H.node[out_node]['parent'] = in_node

        added_nodes += 1
        in_nodes.add(out_node)
        out_nodes.remove(out_node)
        assert added_nodes == len(in_nodes)
        assert len(in_nodes) + len(out_nodes) == G.number_of_nodes()
        for n in in_nodes:
            assert G.has_node(n)
            assert H.has_node(n)

    for u in G.nodes():
        assert H.has_node(u)

    return H
   
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

def centroid_mst(G):
    cent_mst = G.copy()
    cent_mst.remove_edges_from(G.edges())
    
    centroidp = centroid(G)
    cent_mst.add_node('centroid')
    cent_mst.node['centroid']['label'] = 'centroid'
    cent_mst.node['centroid']['coord'] = centroidp
    for u in G.nodes_iter():
        cent_mst.add_edge(u, 'centroid')
        cent_mst[u]['centroid']['length'] = point_dist(cent_mst.node[u]['coord'], centroidp)
    return cent_mst

def main():
    alpha = float(argv[1])
    G = random_point_graph(20)
    mst = pareto_prim(G, alpha)
    steiner = pareto_steiner(G, alpha)
    print graph_costs(mst)
    print graph_costs(steiner)

if __name__ == '__main__':
    main()
