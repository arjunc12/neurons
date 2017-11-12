import networkx
import matplotlib as mpl
mpl.use('agg')
import pylab
from pareto_functions import *
from time import time
from cost_functions import graph_costs, normalize_cost, partially_dominates, mst_cost, satellite_cost
from graph_utils import complete_graph, is_tree
from neuron_utils import pareto_cost, sort_neighbors
from random_graphs import random_point_graph
import argparse
import pandas as pd
from scipy.stats import binom_test

MIN_POINTS = 8
MAX_POINTS = 50

COLNAMES = ['num_points', 'steiner_runtime', 'prim_runtime', 'khuller_runtime',\
            'genetic_runtime', 'prim_comparisons', 'prim_dominated',\
            'khuller_comparisons', 'khuller_dominated', 'genetic_comparisons',\
            'genetic_dominated']

def runtimes_stats():
    df = pd.read_csv('test_runtimes.csv', names=COLNAMES)
    #print "greedy success rate", float(sum(df['dominates'])) / float(sum(df['comparisons']))
    df2 = df.groupby('num_points', as_index=False).agg(pylab.mean)
    pylab.figure()
    points = df2['num_points']
    time1 = df2['steiner_runtime']
    time2 = df2['prim_runtime']
    time3 = df2['khuller_runtime']
    time4 = df2['genetic_runtime']
    pylab.plot(points, time1, c='b', label='steiner')
    pylab.plot(points, time2, c='r', label='prim')
    pylab.plot(points, time3, c='k', label='khuller')
    pylab.plot(points, time4, c='g', label='genetic')
    #pylab.ylim(ymin=-1)
    pylab.legend(loc=2)
    pylab.xlabel('number of points')
    pylab.ylabel('rumtime (minutes)')
    pylab.savefig('test_runtimes/runtimes.pdf', format='pdf')
    pylab.close()

    for algorithm in ['prim', 'khuller', 'genetic']:
        print algorithm
        comparisons = df[algorithm + '_comparisons'].sum()
        dominated = df[algorithm + '_dominated'].sum()
        print float(dominated) / float(comparisons)
        print binom_test(dominated, comparisons)

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
            span_tree = nx.minimum_spanning_tree(G, weight='length')
                   
            opt_scost = float(satellite_cost(sat_tree))
            normalize_scost = lambda cost : normalize_cost(cost, opt_scost)
            opt_mcost = float(mst_cost(span_tree))
            normalize_mcost = lambda cost : normalize_cost(cost, opt_mcost)
            
            sort_neighbors(G)
            
            mcosts1 = []
            mcosts2 = []
            mcosts3 = []
            mcosts4 = []
            
            scosts1 = []
            scosts2 = []
            scosts3 = []
            scosts4 = []

            times1 = []
            times2 = []
            times3 = []
            times4 = []

            delta = 0.01
            alphas = pylab.arange(delta, 1, delta)
            for alpha in alphas:
                print "alpha", alpha
                mcost1, scost1, time1 = time_function(G, alpha, pareto_steiner)
                mcost2, scost2, time2 = time_function(G, alpha, pareto_prim)
                mcost3, scost3, time3 = time_function(G, alpha, pareto_khuller)

                #mcost1 = normalize_mcost(mcost1)
                #scost1 = normalize_scost(scost1)

                mcosts1.append(mcost1)                
                mcosts2.append(mcost2)
                mcosts3.append(mcost3)

                scosts1.append(scost1)
                scosts2.append(scost2)
                scosts3.append(scost3)

                times1.append(time1)
                times2.append(time2)
                times3.append(time3)
            
            steiner_runtime = sum(times1)
            prim_runtime = sum(times2)
            khuller_runtime = sum(times3)

            genetic_start = time()
            genetic_trees = pareto_genetic(G)
            genetic_end = time()

            genetic_runtime = (genetic_end - genetic_start) / 60.0

            for mcost4, scost4, T in genetic_trees: 
                mcosts4.append(mcost4)
                scosts4.append(scost4)

            prim_comparisons, prim_dominated = prop_dominated(mcosts1,\
                                                              scosts1,\
                                                              mcosts2,\
                                                              scosts2)
            khuller_comparisons, khuller_dominated = prop_dominated(mcosts1,\
                                                                    scosts1,\
                                                                    mcosts3,\
                                                                    scosts3)
            genetic_comparisons, genetic_dominated = prop_dominated(mcosts1,\
                                                                    scosts1,\
                                                                    mcosts4,\
                                                                    scosts4)

            pylab.figure()
            pylab.scatter(mcosts1, scosts1, c='b', label='steiner')
            pylab.plot(mcosts1, scosts1, c='b')

            pylab.scatter(mcosts2, scosts2, c='r', label='prim')
            pylab.plot(mcosts2, scosts2, c='r')

            pylab.scatter(mcosts3, scosts3, c='k', label='khuller')
            pylab.plot(mcosts3, scosts3, c='k')

            pylab.scatter(mcosts4, scosts4, c='g', label='genetic')

            pylab.legend()
            pylab.xlabel('spanning tree cost')
            pylab.ylabel('satellite cost')

            pylab.savefig('test_runtimes/pareto_front%d.pdf' % num_points)
            
            pylab.close()

            write_items = [num_points, steiner_runtime, prim_runtime,\
                           khuller_runtime, genetic_runtime, prim_comparisons,\
                           prim_dominated, khuller_comparisons,\
                           khuller_dominated, genetic_comparisons,\
                           genetic_dominated]
            
            write_items = map(str, write_items)
            write_items = ', '.join(write_items)

            outfile = open('test_runtimes.csv', 'a')
            outfile.write('%s\n' % write_items)
            outfile.close()

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
