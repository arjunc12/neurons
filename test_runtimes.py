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

GENETIC = True

def runtimes_stats():
    df = pd.read_csv('test_runtimes.csv', skipinitialspace=True)
    print "total trials"
    print len(df['algorithm']) / len(df['algorithm'].unique())

    ratios = []
    labels = []
    weights = []
    hist_algorithms = ['prim', 'khuller']
    algorithm_labels = {'prim' : 'spanning', 'khuller' : 'khuller'}

    pylab.figure()
    for algorithm, group in df.groupby('algorithm'):
        print algorithm
        comparisons = group['comparisons'].sum()
        dominated = group['dominated'].sum()
        print float(dominated) / float(comparisons)
        print binom_test(dominated, comparisons)
        group = group.groupby('points', as_index=False).agg(pylab.mean)
        pylab.plot(group['points'], group['runtime'], label=algorithm)
        ratio = group['cost ratio']
        ratio = ratio[~pylab.isnan(ratio)]
        print "cost comparisons", len(ratio)
        print "cost ratio", pylab.nanmean(ratio), "+/-", pylab.nanstd(ratio, ddof=1)

        if algorithm in hist_algorithms:
            ratios.append(ratio)
            labels.append(algorithm_labels[algorithm])
            weight = pylab.ones_like(ratio) / float(len(ratio))
            weights.append(weight)

    pylab.legend(loc=2)
    pylab.xlabel('number of points')
    pylab.ylabel('rumtime (minutes)')
    pylab.savefig('test_runtimes/runtimes.pdf', format='pdf')
    pylab.close()

    pylab.figure()
    pylab.hist(ratios, label=labels, weights=weights)
    pylab.xlabel('cost ratio', size=20)
    pylab.ylabel('proportion', size=20)
    pylab.legend()
    pylab.savefig('test_runtimes/cost_ratio_hist.pdf', format='pdf')
    pylab.close()
    
def time_function(G, alpha, pareto_func):
    start = time()
    mst = pareto_func(G, alpha)
    end = time()
    runtime = (end - start) / 60.0
    mcost = mst_cost(mst)
    scost = satellite_cost(mst, relevant_nodes=G.nodes())
    return mcost, scost, runtime

def test_runtimes(num_iters=10, min_points=MIN_POINTS, max_points=MAX_POINTS):    
    for i in xrange(num_iters):
        print "iteration", i
        for num_points in xrange(min_points, max_points + 1):
            print "num_points", num_points
            G = random_point_graph(num_points=num_points)
            
            sat_tree = satellite_tree(G)
            
            sort_neighbors(G)
            
            mcosts = defaultdict(list)
            scosts = defaultdict(list)
            costs = defaultdict(list)
            times = defaultdict(float)

            delta = 0.01
            alphas = pylab.arange(delta, 1, delta)
            algorithms = [pareto_steiner_space, pareto_steiner_space2, pareto_steiner_fast, pareto_steiner_old]
            names = ['space efficient', 'medium space efficient', 'fast', 'unoptimized']
            algorithms += [pareto_prim, pareto_khuller]
            names += ['prim', 'khuller']
            baseline = 'fast'
            for alpha in alphas:
                print "alpha", alpha
                for algorithm, name in zip(algorithms, names):
                    mcost, scost, runtime = time_function(G, alpha, algorithm)
                    cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)
                    mcosts[name].append(mcost)
                    scosts[name].append(scost)
                    costs[name].append(cost)
                    times[name] += runtime

            if num_points <= 50 and GENETIC:
                algorithms.append(pareto_genetic)
                names.append('genetic')
                genetic_start = time()
                genetic_trees = pareto_genetic(G)
                genetic_end = time()

                genetic_runtime = (genetic_end - genetic_start) / 60.0
                times['genetic'] = genetic_runtime
               
                for mcost, scost, T in genetic_trees: 
                    mcosts['genetic'].append(mcost)
                    scosts['genetic'].append(scost)
            
            pylab.figure()
            for name in names:
                mcost = mcosts[name]
                scost = scosts[name]
                pylab.scatter(mcost, scost, label=name)
                pylab.plot(mcost, scost)

            pylab.legend()
            pylab.xlabel('spanning tree cost')
            pylab.ylabel('satellite cost')

            pylab.savefig('test_runtimes/pareto_front%d.pdf' % num_points)
            
            pylab.close()

            header_line = None
            if not os.path.exists('test_runtimes.csv'):
                header_line = ['algorithm', 'points', 'runtime', 'comparisons', 'dominated']
                header_line = ', '.join(header_line)

            mcosts1, scosts1, costs1 = mcosts[baseline], scosts[baseline], pylab.array(costs[baseline])
            with open('test_runtimes.csv', 'a') as outfile:
                if header_line != None:
                    outfile.write('%s\n' % header_line)
                for name in names:
                    write_items = [name, num_points, times[name]]
                    if name in mcosts:
                        assert name in scosts
                        mcosts2, scosts2 = mcosts[name], scosts[name]
                        comparisons, dominated = prop_dominated(mcosts1, scosts1,\
                                                                mcosts2, scosts2)
                        write_items += [comparisons, dominated]
                    else:
                        write_items += ['', '']
                       
                    if name in costs:
                        costs2 = pylab.array(costs[name])
                        cost_ratio = pylab.mean(costs2 / costs1)
                        write_items.append(cost_ratio)
                    else:
                        write_items.append('')
     
                    write_items = map(str, write_items)
                    write_items = ', '.join(write_items)
                    outfile.write('%s\n' % write_items)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', action='store_true')
    parser.add_argument('-n', '--num_iters', type=int, default=10)
    parser.add_argument('-s', '--stats', action='store_true')
    parser.add_argument('--min_points', type=int, default=8)
    parser.add_argument('--max_points', type=int, default=50)

    args = parser.parse_args()
    test = args.test
    num_iters = args.num_iters
    stats = args.stats
    min_points = args.min_points
    max_points = args.max_points
    if test:
        test_runtimes(num_iters=num_iters,\
                      min_points=min_points, max_points=max_points)
    if stats:
        runtimes_stats()

if __name__ == '__main__':
    main()
