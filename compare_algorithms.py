import networkx
import matplotlib as mpl
mpl.use('agg')
import pylab
from pareto_functions import *
from time import time
from cost_functions import graph_costs, normalize_cost, partially_dominates, mst_cost, satellite_cost, pareto_cost
from graph_utils import complete_graph, is_tree
from neuron_utils import sort_neighbors
from random_graphs import random_point_graph_uniform, random_point_graph_normal
import argparse
import pandas as pd
from scipy.stats import binom_test
from collections import defaultdict
import os
import seaborn as sns

MIN_POINTS = 8
MAX_POINTS = 50

BASELINE = 'prim'

def compare_algorithms_stats(model='uniform'):
    df = pd.read_csv('compare_algorithms_%s.csv' % model, skipinitialspace=True)
    algorithms = []
    algorithm_labels = {'khuller' : 'Khuller', 'prim' : 'Karger', 'steiner' : 'Greedy'}
    
    ratios = []
    weights = []
    labels = []
    
    sns.set()
    
    pylab.figure()
    for algorithm, algorithm_group in df.groupby('algorithm'):
        print algorithm
        ratio = algorithm_group['cost ratio']
        print len(ratio) / 99, "point sets", len(ratio), "trials"
        print "success rate", pylab.mean(ratio < 0)
        print "cost ratio", pylab.mean(ratio), "+/-", pylab.std(ratio, ddof=1)
        print "dominated", pylab.mean(algorithm_group['dominated'])

        if algorithm != BASELINE:
            ratios.append(ratio)
            weight = pylab.ones_like(ratio) / float(len(ratio))
            weights.append(weight)
            label = algorithm_labels[algorithm]
            labels.append(label)

        '''
        pylab.subplot(311)
        x = pylab.arange(len(ratio))
        pylab.scatter(x, ratio, label=label)
        pylab.xlabel('test set')
        pylab.ylabel('cost ratio')
        '''

        pylab.subplot(211)
        points = algorithm_group['points']
        pylab.scatter(points, ratio, label='label')
        pylab.xlabel('points')
        pylab.ylabel('cost ratio')

        pylab.subplot(212)
        alpha = algorithm_group['alpha']
        pylab.scatter(alpha, ratio, label=label)
        pylab.xlabel('alpha')
        pylab.ylabel('cost ratio')

    pylab.savefig('compare_algorithms/compare_algorithms_scatter.pdf', format='pdf')
    pylab.close()

    pylab.figure() 
    pylab.hist(ratios, weights=weights, label=labels)
    pylab.legend()
    ax = pylab.gca()
    pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
    pylab.xlabel('Pareto cost (percent higher/lower than greedy)', size=20)
    pylab.ylabel('proportion', size=20)
    pylab.tight_layout()
    pylab.savefig('compare_algorithms/compare_algorithms_hist.pdf', format='pdf')
    pylab.close()

def percent_delta(x, y):
    return 100 * (x - y) / y

def compare_algorithms(num_iters=10, min_points=MIN_POINTS, max_points=MAX_POINTS,\
                       model='uniform'):
    for i in xrange(num_iters):
        print "iteration", i
        for num_points in xrange(min_points, max_points + 1):
            print "num_points", num_points
            G = None
            if model == 'uniform':
                G = random_point_graph_uniform(num_points=num_points)
            elif model == 'normal':
                G = random_point_graph_normal(num_points=num_points)
            else:
                raise ValueError('invalid random graph model')
            
            sort_neighbors(G)
            
            mcosts = defaultdict(list)
            scosts = defaultdict(list)
            costs = defaultdict(list)

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

            baseline_mcosts = mcosts[BASELINE]
            baseline_scosts = scosts[BASELINE]
            baseline_costs = costs[BASELINE]

            header_line = None
            outfname = 'compare_algorithms_%s.csv' % model
            if not os.path.exists(outfname):
                header_line = ['algorithm', 'points', 'alpha', 'mcost', 'scost', 'cost', 'mcost ratio', 'scost ratio', 'cost ratio', 'dominated']
                header_line = ', '.join(header_line)

            with open(outfname, 'a') as outfile:
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

                        baseline_mcost = baseline_mcosts[i]
                        baseline_scost = baseline_scosts[i]
                        baseline_cost = baseline_costs[i]

                        mcost_ratio = percent_delta(mcost, baseline_mcost)
                        scost_ratio = percent_delta(scost, baseline_scost)
                        cost_ratio = percent_delta(cost, baseline_cost)
                        write_items += [mcost_ratio, scost_ratio, cost_ratio]

                        is_dominated = 0
                        for name2 in names:
                            mcosts2 = mcosts[name2]
                            scosts2 = scosts[name2]
                            comparisons, dominated = prop_dominated(mcosts2, scosts2, [mcost], [scost])
                            if dominated > 0:
                                is_dominated = 1
                        write_items.append(is_dominated)

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
    parser.add_argument('--model', choices=['uniform', 'normal'], default='uniform')

    args = parser.parse_args()
    compare = args.compare
    num_iters = args.num_iters
    stats = args.stats
    min_points = args.min_points
    max_points = args.max_points
    model = args.model
    
    if compare:
        compare_algorithms(num_iters=num_iters, min_points=min_points,\
                           max_points=max_points, model=model)
    if stats:
        compare_algorithms_stats(model=model)

if __name__ == '__main__':
    main()
