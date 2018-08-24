import networkx
import matplotlib as mpl
mpl.use('agg')
import pylab
from pareto_functions import *
from time import time
from cost_functions import graph_costs, normalize_cost, partially_dominates, mst_cost, satellite_cost, pareto_cost
from graph_utils import complete_graph, is_tree
from neuron_utils import sort_neighbors
from random_graphs import random_point_graph
import argparse
import pandas as pd
from scipy.stats import binom_test
from collections import defaultdict
import os

MIN_POINTS = 8
MAX_POINTS = 50

def compare_algorithms_stats():
    df = pd.read_csv('compare_algorithms.csv', skipinitialspace=True)
    algorithm_costs = {}

    for name, group in df.groupby('algorithm'):
        algorithm_costs[name] = pylab.array(group['cost'])

    algorithm_labels = {'khuller' : 'khuller', 'prim' : 'Karger', 'steiner' : 'greedy steiner'}

    baseline = algorithm_costs['steiner']
    order = pylab.argsort(baseline)
    baseline = baseline[order]
    x = range(len(order))

    pylab.figure()
    for algorithm in sorted(algorithm_costs.keys()):
        costs = algorithm_costs[algorithm]
        costs = costs[order]
        pylab.scatter(x, 100 * (costs - baseline) / baseline, label=algorithm_labels[algorithm])
        
        print algorithm
        print ((costs < baseline).sum()) / float(len(baseline))

    pylab.legend()
    pylab.xlabel('test set')
    pylab.ylabel('pareto cost vs steiner pareto cost (percent increase/decrease)')
    pylab.savefig('compare_algorithms/compare_algorithms.pdf', format='pdf')
    pylab.close()

def compare_algorithms(num_iters=10, min_points=MIN_POINTS, max_points=MAX_POINTS):    
    for i in xrange(num_iters):
        print "iteration", i
        for num_points in xrange(min_points, max_points + 1):
            print "num_points", num_points
            G = random_point_graph(num_points=num_points)            
            
            sort_neighbors(G)
            
            mcosts = defaultdict(list)
            scosts = defaultdict(list)
            costs = defaultdict(list)
            times = defaultdict(float)

            delta = 0.01
            alphas = pylab.arange(delta, 1, delta)
            algorithms = [pareto_steiner, pareto_prim, pareto_khuller]
            names = ['steiner', 'prim', 'khuller']
            for alpha in alphas:
                print "alpha", alpha
                for algorithm, name in zip(algorithms, names):
                    pareto_tree = algorithm(G, alpha)
                    assert is_tree(pareto_tree)
                    mcost, scost, = graph_costs(pareto_tree, relevant_nodes=G.nodes())
                    cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
                    mcosts[name].append(mcost)
                    scosts[name].append(scost)
                    costs[name].append(cost)
 
            pylab.figure()
            for name in names:
                mcost = mcosts[name]
                scost = scosts[name]
                pylab.scatter(mcost, scost, label=name)
                pylab.plot(mcost, scost)

            pylab.legend()
            pylab.xlabel('spanning tree cost')
            pylab.ylabel('satellite cost')

            pylab.savefig('compare_algorithms/pareto_front%d.pdf' % num_points)
            
            pylab.close()

            header_line = None
            if not os.path.exists('compare_algorithms.csv'):
                header_line = ['algorithm', 'points', 'alpha', 'mcost', 'scost', 'cost']
                header_line = ', '.join(header_line)

            with open('compare_algorithms.csv', 'a') as outfile:
                if header_line != None:
                    outfile.write('%s\n' % header_line)
                for i in xrange(len(alphas)):
                    alpha = alphas[i]
                    for name in names:
                        write_items = [name, num_points, alpha]
                        mcost = mcosts[name][i]
                        scost = scosts[name][i]
                        cost = costs[name][i]
                        write_items += [mcost, scost, cost]

                        write_items = map(str, write_items)
                        write_items = ', '.join(write_items)
                        outfile.write('%s\n' % write_items)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--compare', action='store_true')
    parser.add_argument('-n', '--num_iters', type=int, default=10)
    parser.add_argument('-s', '--stats', action='store_true')
    parser.add_argument('--min_points', type=int, default=8)
    parser.add_argument('--max_points', type=int, default=50)

    args = parser.parse_args()
    compare = args.compare
    num_iters = args.num_iters
    stats = args.stats
    min_points = args.min_points
    max_points = args.max_points
    
    if compare:
        compare_algorithms(num_iters=num_iters,\
                           min_points=min_points, max_points=max_points)
    if stats:
        compare_algorithms_stats()

if __name__ == '__main__':
    main()
