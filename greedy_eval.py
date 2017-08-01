from pareto_functions import pareto_prim, pareto_brute_force
import networkx as nx
from random_graphs import random_point_graph
from neuron_utils import sort_neighbors, pareto_cost
from enumerate_trees import find_all_spanning_trees
import pylab
from cost_functions import graph_costs
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

 
    mcosts = []
    scosts = []
    opt_mcosts = []
    opt_scosts = []
    
    attempts = 0
    successes = 0
    errors = 0
    
    delta = 0.01
    alphas = pylab.arange(0.01, 1, delta)
    for alpha in alphas:
        print "alpha", alpha
        greedy_tree = pareto_prim(G, alpha)
        mcost, scost = graph_costs(greedy_tree)
        mcosts.append(mcost)
        scosts.append(scost)

        opt_tree = pareto_brute_force(G, alpha, trees=span_trees)
        opt_mcost, opt_scost = graph_costs(opt_tree)
        opt_mcosts.append(opt_mcost)
        opt_scosts.append(opt_scost)

        cost = pareto_cost(mcost, scost, alpha)
        opt_cost = pareto_cost(opt_mcost, opt_scost, alpha)

        attempts += 1
        if cost == opt_cost:
            successes += 1
        elif cost < opt_cost:
            errors += 1

    print "attempts", attempts
    print "successes", successes
    print "errors", errors

    num_points = G.number_of_nodes()
    write_items = [num_points, attempts, successes, errors]
    write_items = map(str, write_items)
    write_items = ', '.join(write_items)
    f = open('greedy_eval.csv', 'a')
    f.write('%s\n' % write_items)
    f.close()

    pylab.figure()
    pylab.scatter(mcosts, scosts, c='b', label='greedy')
    pylab.scatter(opt_mcosts, opt_scosts, c='r', label='brute force')

    pylab.plot(mcosts, scosts, c='b')
    pylab.plot(opt_mcosts, opt_scosts, c='r')

    pylab.savefig('greedy_eval/greedy_eval%d.pdf' % num_points)
    pylab.close()

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
