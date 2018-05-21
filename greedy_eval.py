from pareto_functions import *
import networkx as nx
from random_graphs import random_point_graph
from neuron_utils import *
from enumerate_trees import find_all_spanning_trees
import pylab
from cost_functions import *
from sys import argv
import argparse
import pandas as pd
from collections import defaultdict
import os

def greedy_eval(G):
    print "sorting neighbors"
    sort_neighbors(G)
    print "getting trees"
    all_trees = find_all_spanning_trees(G, root=G.graph['root'])
    
    span_trees = []
    for tree in all_trees:
        span_tree = G.copy()
        span_tree.remove_edges_from(G.edges())
        for u, v in tree.edges():
            span_tree.add_edge(u, v)
            span_tree[u][v] = G[u][v]
            span_tree[v][u] = G[v][u]
        span_trees.append(span_tree)

    mst_cost = best_mst_cost(G)
    sat_cost = best_satellite_cost(G)
 
    mcosts = []
    scosts = []
    opt_mcosts = []
    opt_scosts = []
    
    attempts = 0
    successes = 0
    errors = 0
    
    delta = 0.01
    alphas = pylab.arange(0.01, 1, delta)
    fname = 'greedy_eval.csv'
    first_time = not os.path.exists(fname)
    f = open('greedy_eval.csv', 'a')
    if first_time:
        f.write('algorithm, points, alpha, mcost, scost, cost, optimal mcost, ' +\
                 'optimal scost, optimal cost, success, error\n')
    num_points = G.number_of_nodes()
    
    prev_mcost = None
    prev_scost = None
    
    prev_opt_mcost = None
    prev_opt_scost = None
    
    algorithms = [pareto_prim, pareto_khuller, pareto_steiner]
    names = ['prim', 'khuller', 'steiner']
    
    for alpha in alphas:
        print "alpha", alpha
        attempts += 1

        opt_tree = pareto_brute_force(G, alpha, trees=span_trees)
        opt_mcost, opt_scost = graph_costs(opt_tree, relevant_nodes=G.nodes())
        opt_cost = pareto_cost(opt_mcost, opt_scost, alpha)
        
        for algorithm, name in zip(algorithms, names):
            tree = algorithm(G, alpha)
            mcost, scost = graph_costs(tree, relevant_nodes=G.nodes())
            cost = pareto_cost(mcost, scost, alpha)

            success = cost == opt_cost
            error = cost < opt_cost
            if error:
                pass
                '''
                print "error", algorithm
                print "mcost", mcost, opt_mcost
                print "scost", scost, opt_scost
                print "cost", cost, opt_cost
                '''
    
            write_items = [name, num_points, alpha, mcost, scost, cost,\
                           opt_mcost, opt_scost, opt_cost,\
                           int(success), int(error)]
            write_items = map(str, write_items)
            write_items = ', '.join(write_items)
            f.write('%s\n' % write_items)

    f.close()

def greedy_eval_stats():
    df = pd.read_csv('greedy_eval.csv', skipinitialspace=True)

    for algorithm, group in df.groupby('algorithm'):
        print algorithm 
        print "success rate", pylab.mean(group['success'])
        print "error rate", pylab.mean(group['error'])
        #print group[group['error'] == 1]
        
        print "mcost ratio", pylab.mean(group['mcost'] / group['optimal mcost'])
        print "scost ratio", pylab.mean(group['scost'] / group['optimal scost'])
        print "cost ratio", pylab.mean(group['cost'] / group['optimal cost'])
        
        print "max mcost ratio", max(group['mcost'] / group['optimal mcost'])
        print "max scost ratio", max(group['scost'] / group['optimal scost'])
        print "max cost ratio", max(group['cost'] / group['optimal cost'])

    pylab.figure()
    pylab.scatter(df['mcost'], df['optimal mcost'])
    pylab.savefig('greedy_eval/greedy_eval_mcost.pdf', format='pdf')
    pylab.close()

    pylab.figure()
    pylab.scatter(df['scost'], df['optimal scost'])
    pylab.savefig('greedy_eval/greedy_eval_scost.pdf', format='pdf')
    pylab.close()

    pylab.figure()
    pylab.scatter(df['cost'], df['optimal cost'])
    pylab.savefig('greedy_eval/greedy_eval_cost.pdf', format='pdf')
    pylab.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--min_points', type=int, default=5)
    parser.add_argument('--max_points', type=int, default=9)
    parser.add_argument('-x', '--num_iters', type=int, default=10)
    parser.add_argument('-e', '--evaluate', action='store_true')
    parser.add_argument('-s', '--stats', action='store_true')

    args = parser.parse_args()
    min_points = args.min_points
    max_points = args.max_points
    num_iters = args.num_iters
    evaluate = args.evaluate
    stats = args.stats

    if stats:
        greedy_eval_stats()

    if evaluate:
        for i in xrange(num_iters):
            for num_points in xrange(min_points, max_points + 1):
                G = random_point_graph(num_points=num_points)
                greedy_eval(G)

if __name__ == '__main__':
    main()
