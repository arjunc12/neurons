import networkx as nx
from neuron_utils import *
from kruskal import kruskal
from random_graphs import random_mst
from itertools import combinations
import numpy as np
from cost_functions import *
import pylab
from bisect import bisect_left, insort
from random import sample, choice, uniform, randint, random, seed
from collections import defaultdict
from prufer import *
from enumerate_trees import find_all_spanning_trees
from random_graphs import random_point_graph
from steiner_midpoint import *
from sys import argv
from scipy.spatial.distance import euclidean
from graph_utils import *
from khuller import khuller
from time import time
import argparse

POP_SIZE = 400
GENERATIONS = 20000
MUTATION_PROB = 0.01

STEINER_MIDPOINTS = 10

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

def pareto_steiner(G, alpha, axon=False):
    return pareto_steiner_fast(G, alpha, axon=axon)

def pareto_steiner_space(G, alpha, axon=False):
    root = G.graph['root']

    H = nx.Graph()
   
    H.add_node(root)
    H.graph['root'] = root
    H.node[root]['droot'] = 0
    H.node[root]['parent'] = None
    root_coord = G.node[root]['coord']
    H.node[root]['coord'] = root_coord
    H.node[root]['label'] = 'root'
    added_nodes = 1

    out_nodes = set(G.nodes())
    out_nodes.remove(root)
   
    graph_mcost = 0
    graph_scost = 0

    unpaired_nodes = [root]

    node_index = max(G.nodes()) + 1

    dist_error = 0

    steps = 0

    best_edges = []
    while added_nodes < G.number_of_nodes():
        for u in unpaired_nodes:
            if axon and (u == H.graph['root']) and (H.degree(u) > 0):
                continue
        
            closest_neighbor = None
            closest_dist = float("inf")
            for v in out_nodes:
                coord1 = H.node[u]['coord']
                coord2 = G.node[v]['coord']
                dist = point_dist(coord1, coord2)
                if dist < closest_dist:
                    closest_dist = dist
                    closest_neighbor = v
 
            p1 = H.node[u]['coord']
            p2 = G.node[closest_neighbor]['coord']
                
            length = point_dist(p1, p2)
            mcost = length
            scost = length + H.node[u]['droot']
            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
            insort(best_edges, (cost, u, closest_neighbor))

        cost, u, v = best_edges.pop(0)

        best_edges2 = []
        unpaired_nodes = [u, v]
        while len(best_edges) > 0:
            cost, x, y = best_edges.pop(0)
            if y == v:
                unpaired_nodes.append(x)
            else:
                best_edges2.append((cost, x, y))
        best_edges = best_edges2

        H.add_node(v)
        H.node[v]['coord'] = G.node[v]['coord']
        H.node[v]['label'] = 'synapse'
        out_nodes.remove(v)

        p1 = H.node[u]['coord']
        p2 = H.node[v]['coord']
        midpoints = steiner_points(p1, p2, npoints=STEINER_MIDPOINTS)
        midpoint_nodes = []
        for midpoint in midpoints:
            midpoint_node = node_index
            node_index += 1
            H.add_node(midpoint_node)
            H.node[midpoint_node]['coord'] = midpoint

            midpoint_nodes.append(midpoint_node)

            unpaired_nodes.append(midpoint_node)

        line_nodes = [v] + list(reversed(midpoint_nodes)) + [u]
        for i in xrange(-1, -len(line_nodes), -1):
            n1 = line_nodes[i]
            n2 = line_nodes[i - 1]
            H.add_edge(n1, n2)
            H[n1][n2]['length'] = node_dist(H, n1, n2)
            H.node[n2]['parent'] = n1
            H.node[n2]['droot'] = node_dist(H, n2, u) + H.node[u]['droot']
            if not G.has_node(n2):
                H.node[n2]['label'] = 'steiner_midpoint'

        added_nodes += 1

    return H

def pareto_steiner_space2(G, alpha, axon=False):
    root = G.graph['root']

    H = nx.Graph()
   
    H.add_node(root)
    H.graph['root'] = root
    H.node[root]['droot'] = 0
    H.node[root]['parent'] = None
    root_coord = G.node[root]['coord']
    H.node[root]['coord'] = root_coord
    H.node[root]['label'] = 'root'
    added_nodes = 1

    out_nodes = set(G.nodes())
    out_nodes.remove(root)
   
    graph_mcost = 0
    graph_scost = 0

    closest_neighbors = {}
    is_sorted = 'sorted' in G.graph
    max_neighbors = 10
    for u in G.nodes_iter():
        if is_sorted:
            closest_neighbors[u] = G.node[u]['close_neighbors'][:max_neighbors]
        else:
            closest_neighbors[u] = k_nearest_neighbors(G, u, k=max_neighbors, candidate_nodes=None)
    best_neighbor = {}

    unpaired_nodes = set([root])

    node_index = max(G.nodes()) + 1

    dist_error = 0

    steps = 0

    best_edges = []
    while added_nodes < G.number_of_nodes():
        assert len(out_nodes) > 0
        best_edge = None
        best_mcost = None
        best_scost = None
        best_cost = float("inf")

        best_choice = None
        best_midpoint = None

        candidate_edges = []
        for u in unpaired_nodes:
            if axon and (u == H.graph['root']) and (H.degree(u) > 0):
                continue

            assert H.has_node(u)
            assert 'droot' in H.node[u]
        
            closest_neighbor = None
            while closest_neighbor == None:
                invalid_neighbors = []
            
                for i in xrange(len(closest_neighbors[u])):
                    v = closest_neighbors[u][i]
                    if H.has_node(v):
                        invalid_neighbors.append(v)
                    else:
                        closest_neighbor = v
                        break

                for invalid_neighbor in invalid_neighbors:
                    closest_neighbors[u].remove(invalid_neighbor)

                if closest_neighbor == None:
                    closest_neighbors[u] = k_nearest_neighbors(G, u, k=10, candidate_nodes=out_nodes)
            
            
            assert closest_neighbor != None
            assert not H.has_node(closest_neighbor)
            
            p1 = H.node[u]['coord']
            p2 = G.node[closest_neighbor]['coord']
                
            length = point_dist(p1, p2)
            mcost = length
            scost = length + H.node[u]['droot']
            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
            insort(best_edges, (cost, u, closest_neighbor))

        cost, u, v = best_edges.pop(0)

        best_edges2 = []
        unpaired_nodes = set([u, v])
        for cost, x, y in best_edges:
            if y == v:
                unpaired_nodes.add(x)
            else:
                best_edges2.append((cost, x, y))
        best_edges = best_edges2

        assert H.has_node(u)
        assert not H.has_node(v)
        H.add_node(v)
        H.node[v]['coord'] = G.node[v]['coord']
        H.node[v]['label'] = 'synapse'
        out_nodes.remove(v)

        p1 = H.node[u]['coord']
        p2 = H.node[v]['coord']
        midpoints = steiner_points(p1, p2, npoints=STEINER_MIDPOINTS)
        midpoint_nodes = []
        for midpoint in midpoints:
            midpoint_node = node_index
            node_index += 1
            H.add_node(midpoint_node)
            H.node[midpoint_node]['coord'] = midpoint

            neighbors = []
            for out_node in out_nodes:
                out_coord = G.node[out_node]['coord']
                dist = point_dist(midpoint, out_coord)
                neighbors.append((dist, out_node))

            neighbors = sorted(neighbors)
            closest_neighbors[midpoint_node] = []
            for dist, neighbor in neighbors:
                closest_neighbors[midpoint_node].append(neighbor)

            midpoint_nodes.append(midpoint_node)

            unpaired_nodes.add(midpoint_node)

        line_nodes = [v] + list(reversed(midpoint_nodes)) + [u]
        for i in xrange(-1, -len(line_nodes), -1):
            n1 = line_nodes[i]
            n2 = line_nodes[i - 1]
            H.add_edge(n1, n2)
            H[n1][n2]['length'] = node_dist(H, n1, n2)
            assert 'droot' in H.node[n1]
            H.node[n2]['parent'] = n1
            H.node[n2]['droot'] = node_dist(H, n2, u) + H.node[u]['droot']
            if not G.has_node(n2):
                H.node[n2]['label'] = 'steiner_midpoint'

        added_nodes += 1
    return H

def pareto_steiner_fast(G, alpha, axon=False):
    root = G.graph['root']

    H = nx.Graph()
   
    H.add_node(root)
    H.graph['root'] = root
    H.node[root]['droot'] = 0
    H.node[root]['parent'] = None
    root_coord = G.node[root]['coord']
    H.node[root]['coord'] = root_coord
    H.node[root]['label'] = 'root'
    added_nodes = 1

    in_nodes = set()
    out_nodes = set(G.nodes())
    in_nodes.add(root)
    out_nodes.remove(root)
   
    graph_mcost = 0
    graph_scost = 0

    closest_neighbors = {}
    is_sorted = 'sorted' in G.graph
    for u in G.nodes_iter():
        if is_sorted:
            closest_neighbors[u] = G.node[u]['close_neighbors'][:]
        else:
            closest_neighbors[u] = k_nearest_neighbors(G, u, k=None, candidate_nodes=None)

    unpaired_nodes = set([root])

    node_index = max(G.nodes()) + 1

    dist_error = 0

    steps = 0

    best_edges = []
    while added_nodes < G.number_of_nodes():
        assert len(out_nodes) > 0
        best_edge = None
        best_mcost = None
        best_scost = None
        best_cost = float("inf")

        best_choice = None
        best_midpoint = None

        candidate_edges = []
        for u in unpaired_nodes:
            if axon and (u == H.graph['root']) and (H.degree(u) > 0):
                continue

            assert H.has_node(u)
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

            for invalid_neighbor in invalid_neighbors:
                closest_neighbors[u].remove(invalid_neighbor)
           
            assert closest_neighbor != None
            assert not H.has_node(closest_neighbor)
            
            p1 = H.node[u]['coord']
            p2 = G.node[closest_neighbor]['coord']
                
            length = point_dist(p1, p2)
            mcost = length
            scost = length + H.node[u]['droot']
            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
            insort(best_edges, (cost, u, closest_neighbor))

        cost, u, v = best_edges.pop(0)

        best_edges2 = []
        unpaired_nodes = set([u, v])
        for cost, x, y in best_edges:
            if y == v:
                unpaired_nodes.add(x)
            else:
                best_edges2.append((cost, x, y))
        best_edges = best_edges2

        assert H.has_node(u)
        assert not H.has_node(v)
        H.add_node(v)
        H.node[v]['coord'] = G.node[v]['coord']
        H.node[v]['label'] = 'synapse'
        in_nodes.add(v)
        out_nodes.remove(v)

        p1 = H.node[u]['coord']
        p2 = H.node[v]['coord']
        midpoints = steiner_points(p1, p2, npoints=STEINER_MIDPOINTS)
        midpoint_nodes = []
        for midpoint in midpoints:
            midpoint_node = node_index
            node_index += 1
            H.add_node(midpoint_node)
            H.node[midpoint_node]['coord'] = midpoint

            neighbors = []
            for out_node in out_nodes:
                out_coord = G.node[out_node]['coord']
                dist = point_dist(midpoint, out_coord)
                neighbors.append((dist, out_node))

            neighbors = sorted(neighbors)
            closest_neighbors[midpoint_node] = []
            for dist, neighbor in neighbors:
                closest_neighbors[midpoint_node].append(neighbor)

            midpoint_nodes.append(midpoint_node)

            unpaired_nodes.add(midpoint_node)

        line_nodes = [v] + list(reversed(midpoint_nodes)) + [u]
        for i in xrange(-1, -len(line_nodes), -1):
            n1 = line_nodes[i]
            n2 = line_nodes[i - 1]
            H.add_edge(n1, n2)
            H[n1][n2]['length'] = node_dist(H, n1, n2)
            assert 'droot' in H.node[n1]
            H.node[n2]['parent'] = n1
            H.node[n2]['droot'] = node_dist(H, n2, u) + H.node[u]['droot']
            if not G.has_node(n2):
                H.node[n2]['label'] = 'steiner_midpoint'

        added_nodes += 1
    return H

def pareto_steiner_old(G, alpha, axon=False):
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
    H.node[root]['label'] = 'root'
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
            assert H.has_node(u)
            assert not H.has_node(v)

            p1 = H.node[u]['coord']
            p2 = G.node[v]['coord']
                
            length = point_dist(p1, p2)
            mcost = length
            scost = length + H.node[u]['droot']
            cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
            
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
                best_mcost = mcost
                best_scost = scost 

        if best_edge == None:
            break
        u, v = best_edge
        assert H.has_node(u)
        assert not H.has_node(v)
        H.add_node(v)
        H.node[v]['coord'] = G.node[v]['coord']
        H.node[v]['label'] = 'synapse'
        in_nodes.add(v)
        out_nodes.remove(v)

        p1 = H.node[u]['coord']
        p2 = H.node[v]['coord']
        midpoints = steiner_points(p1, p2, npoints=STEINER_MIDPOINTS)
        midpoint_nodes = []
        for midpoint in midpoints:
            midpoint_node = node_index
            node_index += 1
            H.add_node(midpoint_node)
            H.node[midpoint_node]['coord'] = midpoint

            neighbors = []
            for out_node in out_nodes:
                out_coord = G.node[out_node]['coord']
                dist = point_dist(midpoint, out_coord)
                neighbors.append((dist, out_node))

            neighbors = sorted(neighbors)
            closest_neighbors[midpoint_node] = []
            for dist, neighbor in neighbors:
                closest_neighbors[midpoint_node].append(neighbor)

            midpoint_nodes.append(midpoint_node)

        line_nodes = [v] + list(reversed(midpoint_nodes)) + [u]
        for i in xrange(-1, -len(line_nodes), -1):
            n1 = line_nodes[i]
            n2 = line_nodes[i - 1]
            H.add_edge(n1, n2)
            H[n1][n2]['length'] = node_dist(H, n1, n2)
            assert 'droot' in H.node[n1]
            H.node[n2]['parent'] = n1
            H.node[n2]['droot'] = node_dist(H, n2, u) + H.node[u]['droot']
            if not G.has_node(n2):
                H.node[n2]['label'] = 'steiner_midpoint'

        added_nodes += 1
    return H

def alpha_to_beta(alpha, opt_mcost, opt_scost):
    assert alpha > 0
    assert alpha < 1
    num = 2 * alpha * opt_mcost
    denom = (1 - alpha) * opt_scost
    frac = num / denom
    frac **= 0.5
    return 1 + frac

def pareto_khuller(G, alpha, span_tree=None, sat_tree=None):
    if span_tree == None:
        span_tree = nx.minimum_spanning_tree(G, weight='length')
    if sat_tree == None:
        sat_tree = satellite_tree(G)

    opt_mcost = mst_cost(span_tree)
    opt_scost = satellite_cost(sat_tree)
    beta = alpha_to_beta(alpha, opt_mcost, opt_scost)
    return khuller(G, span_tree, sat_tree, beta)
        

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

def centroid(G):
    root = G.graph['root']
    root_coord = G.node[root]['coord']
    centroid = np.zeros(len(root_coord))
    for u in G.nodes_iter():
        point = G.node[u]['coord']
        assert len(point) == len(root_coord)
        if u != root:
            centroid += point
    centroid /= G.number_of_nodes() - 1
    return centroid

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
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--points', type=int, default=100)
    parser.add_argument('-a', '--alpha', type=float, default=0.5)

    args = parser.parse_args()
    points = args.points
    alpha = args.alpha

    G = random_point_graph(points)

    #algorithms = [pareto_steiner_space, pareto_steiner_space2, pareto_steiner_fast, pareto_steiner_old]
    algorithms = [pareto_steiner_space2, pareto_steiner_fast]
    for pareto_func in algorithms:
        start = time()
        tree = pareto_func(G, alpha)
        end = time()
        #print graph_costs(tree, relevant_nodes=G.nodes())
        print end - start 

if __name__ == '__main__':
   main()
