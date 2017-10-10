from pareto_functions import pareto_prim, pareto_brute_force
import networkx as nx
from random_graphs import random_point_graph
from neuron_utils import sort_neighbors, pareto_cost
from enumerate_trees import find_all_spanning_trees
import pylab
from cost_functions import *
from sys import argv
import argparse
import pandas as pd

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

    normalize_mcost, normalize_scost = graph_normalize_functions(G)

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
    f = open('greedy_eval.csv', 'a')
    num_points = G.number_of_nodes()
    
    prev_mcost = None
    prev_scost = None
    
    prev_opt_mcost = None
    prev_opt_scost = None

    for alpha in alphas:
        print "alpha", alpha
        greedy_tree = pareto_prim(G, alpha)
        mcost, scost = graph_costs(greedy_tree)
        #mcost = normalize_mcost(mcost)
        #scost = normalize_scost(scost)
        mcosts.append(mcost)
        scosts.append(scost)

        opt_tree = pareto_brute_force(G, alpha, trees=span_trees)
        opt_mcost, opt_scost = graph_costs(opt_tree)
        #opt_mcost = normalize_mcost(opt_mcost)
        #opt_scost = normalize_scost(opt_scost)
        opt_mcosts.append(opt_mcost)
        opt_scosts.append(opt_scost)

        cost = pareto_cost(mcost, scost, alpha)
        opt_cost = pareto_cost(opt_mcost, opt_scost, alpha)

        attempts += 1
        if cost == opt_cost:
            successes += 1
        elif cost < opt_cost:
            print "mcost", mcost, opt_mcost, "scost", scost, opt_scost, "cost", cost, opt_cost
            errors += 1
            assert False
    
        write_items = [num_points, alpha, mcost, scost, opt_mcost, opt_scost,\
                      mst_cost, sat_cost]
        write_items = map(str, write_items)
        write_items = ', '.join(write_items)
        f.write('%s\n' % write_items)

    print "attempts", attempts
    print "successes", successes
    print "errors", errors

    #write_items = [num_points, attempts, successes, errors]
    #write_items = map(str, write_items)
    #write_items = ', '.join(write_items)
    #f.write('%s\n' % write_items)
    f.close()

    pylab.figure()
    pylab.scatter(mcosts, scosts, c='b', label='greedy')
    pylab.scatter(opt_mcosts, opt_scosts, c='r', label='brute force')

    pylab.plot(mcosts, scosts, c='b')
    pylab.plot(opt_mcosts, opt_scosts, c='r')

    pylab.savefig('greedy_eval/greedy_eval%d.pdf' % num_points)
    pylab.close()

def greedy_eval_stats():
    df = pd.read_csv('greedy_eval.csv', names=['num_points', 'alpha',\
                                               'mcost', 'scost',\
                                               'opt_mcost', 'opt_scost'])

    print "total trials"
    print len(df['alpha']) / len(df['alpha'].unique())

    df['cost'] = (df['alpha'] * df['mcost']) + ((1 - df['alpha']) * df['scost'])
    df['opt_cost'] = (df['alpha'] * df['opt_mcost']) + ((1 - df['alpha']) * df['opt_scost'])
    
    df['success'] = df['cost'] == df['opt_cost']
    df['success'] = df['success'].astype(int)
    
    df['error'] = df['cost'] < df['opt_cost']
    df['error'] = df['error'].astype(int)
    
    print "errors", sum(df['error'])
    print df[df['error'] == 1]
    
    print "success rate", pylab.mean(df['success'])

    df['mcost_ratio'] = df['mcost'] / df['opt_mcost']
    df['scost_ratio'] = df['scost'] / df['opt_scost']
    df['cost_ratio'] = df['cost'] / df['opt_cost']

    print "mcost ratio", pylab.mean(df['mcost_ratio'])
    print "scost ratio", pylab.mean(df['scost_ratio'])
    print "cost ratio", pylab.mean(df['cost_ratio'])

    print "max mcost_ratio", max(df['mcost_ratio'])
    print "max scost ratio", max(df['scost_ratio'])
    print "max cost ratio", max(df['cost_ratio'])

    pylab.figure()
    pylab.scatter(df['mcost'], df['opt_mcost'])
    pylab.savefig('greedy_eval/greedy_eval_mcost.pdf', format='pdf')
    pylab.close()

    pylab.figure()
    pylab.scatter(df['scost'], df['opt_scost'])
    pylab.savefig('greedy_eval/greedy_eval_scost.pdf', format='pdf')
    pylab.close()

    pylab.figure()
    pylab.scatter(df['cost'], df['opt_cost'])
    pylab.savefig('greedy_eval/greedy_eval_cost.pdf', format='pdf')
    pylab.close()


'''
def greedy_eval_stats():
    df = pd.read_csv('greedy_eval.csv', names = ['num_points', 'attempts', 'successes', 'errors'])
    print "errors", sum(df['errors'])
    df['success_rate'] = df['successes'] / df['attempts']
    pylab.figure()
    pylab.scatter(df['num_points'], df['success_rate'])
    pylab.xlabel('number of points')
    pylab.ylabel('greedy success rate')
    pylab.savefig('greedy_eval/greedy_eval.pdf', format='pdf')
    pylab.close()
'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--min_points', type=int, default=5)
    parser.add_argument('--max_points', type=int, default=10)
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
