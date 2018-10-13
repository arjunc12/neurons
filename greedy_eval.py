from pareto_functions import *
import networkx as nx
from random_graphs import random_point_graph_uniform, random_point_graph_normal
from neuron_utils import *
from enumerate_trees import find_all_spanning_trees
import pylab
from cost_functions import *
from sys import argv
import argparse
import pandas as pd
from collections import defaultdict
import os
import seaborn as sns
from time import time

MIN_POINTS = 5
MAX_POINTS = 8

BASELINE = 'brute force'

def percent_delta(x, y):
    return 100 * (x - y) / y

def greedy_eval(num_iters=10, min_points=MIN_POINTS, max_points=MAX_POINTS, model='uniform'):
    for i in xrange(num_iters):
        for num_points in xrange(min_points, max_points + 1):
            G = None
            if model == 'uniform':
                G = random_point_graph_uniform(num_points=num_points)
            elif model == 'normal':
                G = random_point_graph_normal(num_points=num_points)
            
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
                    for item, val in G[u][v].iteritems():
                        span_tree[u][v][item] = val
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
            alphas = pylab.arange(delta, 1, delta)
            
            prev_mcost = None
            prev_scost = None
            
            prev_opt_mcost = None
            prev_opt_scost = None
            
            def brute_force(G, alpha):
                return pareto_brute_force(G, alpha, trees=span_trees)

            algorithms = [brute_force, pareto_prim, pareto_khuller, pareto_steiner]
            names = ['brute force', 'prim', 'khuller', 'steiner']
            mcosts = defaultdict(list)
            scosts = defaultdict(list)
            costs = defaultdict(list)
            times = defaultdict(list)
            
            for alpha in alphas:
                print "alpha", alpha
                attempts += 1
 
                for algorithm, name in zip(algorithms, names):
                    start = time()
                    tree = algorithm(G, alpha)
                    end = time()
                    mcost, scost = graph_costs(tree, relevant_nodes=G.nodes())
                    cost = pareto_cost(mcost, scost, alpha)
                    algorithm_time = end - start

                    mcosts[name].append(mcost)
                    scosts[name].append(scost)
                    costs[name].append(cost)
                    times[name].append(algorithm_time)

            baseline_mcosts = mcosts[BASELINE]
            baseline_scosts = scosts[BASELINE]
            baseline_costs = costs[BASELINE]
            baseline_times = times[BASELINE]

            header_line = None
            outfname = 'greedy_eval_%s.csv' % model
            if not os.path.exists(outfname):
                header_line = ['algorithm', 'points', 'alpha', 'mcost', 'scost', 'cost', 'time',\
                               'mcost ratio', 'scost ratio', 'cost ratio', 'time ratio', 'dominated']
                header_line = ', '.join(header_line)

            with open(outfname, 'a') as outfile:
                if header_line != None:
                    outfile.write('%s\n' % header_line)
                for i in xrange(len(alphas)):
                    for name in names:
                        alpha = alphas[i]
                        write_items = [name, num_points, alpha]

                        mcost = mcosts[name][i]
                        scost = scosts[name][i]
                        cost = costs[name][i]
                        algorithm_time = times[name][i]
                        write_items += [mcost, scost, cost, algorithm_time]

                        baseline_mcost = baseline_mcosts[i]
                        baseline_scost = baseline_scosts[i]
                        baseline_cost = baseline_costs[i]
                        baseline_time = baseline_times[i]

                        mcost_ratio = percent_delta(mcost, baseline_mcost)
                        scost_ratio = percent_delta(scost, baseline_scost)
                        cost_ratio = percent_delta(cost, baseline_cost)
                        time_ratio = percent_delta(algorithm_time, baseline_time)
                        write_items += [mcost_ratio, scost_ratio, cost_ratio, time_ratio]

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

def greedy_eval_stats(model='uniform'):
    df = pd.read_csv('greedy_eval_%s.csv' % model, skipinitialspace=True)
    df['success'] = df['cost ratio'] == 0
    df['error'] = df['cost ratio'] < 0

    mcost_ratios = []
    scost_ratios = []
    cost_ratios = []
    time_ratios = []
    algorithm_labels = {'brute force' : 'Brute force', 'steiner' : 'Greedy', 'prim' : 'Karger', 'khuller' : 'Khuller'}
    labels = []

    for algorithm, group in df.groupby('algorithm'):
        print algorithm
        alpha = group['alpha']
        print len(alpha) / len(alpha.unique()), "point sets", len(alpha), "trials"
        print "success rate", pylab.mean(group['success'])
        print "error rate", pylab.mean(group['error'])
        print "dominated", pylab.mean(group['dominated'])
        
        mcost_ratio = group['mcost ratio']
        scost_ratio = group['scost ratio']
        cost_ratio = group['cost ratio']
        time_ratio = group['time ratio']

        print "mcost ratio", pylab.nanmean(mcost_ratio), "+/-", pylab.nanstd(mcost_ratio, ddof=1)
        print "scost ratio", pylab.nanmean(scost_ratio), "+/-", pylab.nanstd(scost_ratio, ddof=1)
        print "cost ratio", pylab.nanmean(cost_ratio), "+/-", pylab.nanstd(cost_ratio, ddof=1)
        print "time ratio", pylab.nanmean(time_ratio), "+/-", pylab.nanstd(time_ratio, ddof=1)
        
        print "max mcost ratio", pylab.nanmax(mcost_ratio)
        print "max scost ratio", pylab.nanmax(scost_ratio)
        print "max cost ratio", pylab.nanmax(cost_ratio)
        print "max time ratio", pylab.nanmax(time_ratio)

        mcost_ratios.append(mcost_ratio)
        scost_ratios.append(scost_ratio)
        cost_ratios.append(cost_ratio)
        time_ratios.append(time_ratio)
        labels.append(algorithm_labels[algorithm])

        for points, points_group in group.groupby('points'):
            print points, pylab.mean(points_group['time'])

    sns.set()
    
    axis_text = {'mcost' : 'wiring cost', 'scost' : 'conduction delay', 'cost' : 'Pareto cost', 'time' : 'runtime'}
    for ratios, fname in zip([mcost_ratios, scost_ratios, cost_ratios, time_ratios], ['mcost', 'scost', 'cost', 'time']):
        pylab.figure()
        x = []
        weights = []
        for ratio, label in zip(ratios, labels):
            ratio = pylab.array(ratio)
            ratio = ratio[~pylab.isnan(ratio)]
            x.append(ratio)
            wt = pylab.ones_like(ratio) / float(len(ratio))
            weights.append(wt)
        pylab.hist(x, weights=weights, label=labels)
        pylab.legend()
        ax = pylab.gca()
        pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
        ax_text = axis_text[fname]
        pylab.xlabel('%s (percent higher/lower than %s)' % (ax_text, BASELINE), size=20)
        pylab.ylabel('proportion', size=20)
        pylab.tight_layout()
        pylab.savefig('greedy_eval/%s_ratios_hist.pdf' % fname, format='pdf')
        pylab.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--min_points', type=int, default=MIN_POINTS)
    parser.add_argument('--max_points', type=int, default=MAX_POINTS)
    parser.add_argument('-x', '--num_iters', type=int, default=10)
    parser.add_argument('-e', '--evaluate', action='store_true')
    parser.add_argument('-s', '--stats', action='store_true')
    parser.add_argument('-m', '--model', choices=['uniform', 'normal'], default='uniform')

    args = parser.parse_args()
    min_points = args.min_points
    max_points = args.max_points
    num_iters = args.num_iters
    evaluate = args.evaluate
    stats = args.stats
    model = args.model

    if stats:
        greedy_eval_stats(model)

    if evaluate:
        greedy_eval(num_iters, min_points, max_points, model)

if __name__ == '__main__':
    main()
