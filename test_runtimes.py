import networkx
import matplotlib as mpl
mpl.use('agg')
import pylab
from pareto_functions import *
from time import time
from cost_functions import graph_costs, normalize_cost
from neuron_utils import complete_graph, pareto_cost, sort_neighbors, is_tree
from random_graphs import random_point_graph
import argparse
import pandas as pd

MIN_POINTS = 8
MAX_POINTS = 50

COLNAMES = ['num_points', 'time1', 'time2'] #, 'comparisons', 'dominates']

def runtimes_stats():
    df = pd.read_csv('test_runtimes.csv', names=COLNAMES)
    #print "greedy success rate", float(sum(df['dominates'])) / float(sum(df['comparisons']))
    df = df[['num_points', 'time1', 'time2']]
    df = df.groupby('num_points', as_index=False).agg(pylab.mean)
    pylab.figure()
    time1 = df['time1']
    time2 = df['time2']
    pylab.plot(df['num_points'], time1, c='b', label='greedy')
    pylab.scatter(df['num_points'], time2, c='r', label='genetic')
    #pylab.ylim(ymin=-1)
    pylab.legend(loc=2)
    pylab.xlabel('number of points')
    pylab.ylabel('rumtime (minutes)')
    pylab.savefig('test_runtimes/runtimes.pdf', format='pdf')
    pylab.close()

def time_function(G, alpha, pareto_func):
    start = time()
    mst = pareto_func(G, alpha)
    end = time()
    runtime = (end - start) / 60.0
    mcost, scost = graph_costs(mst)
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
            
            scosts1 = []
            scosts2 = []

            times1 = []
            times2 = []

            delta = 0.01
            comparisons = 0
            dominates = 0
            for alpha in pylab.arange(delta, 1, delta):
                print "alpha", alpha
                mcost1, scost1, time1 = time_function(G, alpha, pareto_prim)

                mcost1 = normalize_mcost(mcost1)
                scost1 = normalize_scost(scost1)

                mcosts1.append(mcost1)                
                scosts1.append(scost1)

                times1.append(time1)

            
            greedy_runtime = sum(times1)

            genetic_start = time()
            genetic_trees = pareto_genetic(G)
            genetic_end = time()

            genetic_runtime = (genetic_end - genetic_start) / 60.0

            for mcost2, scost2, T in genetic_trees:
                mcost2 = normalize_mcost(mcost2)
                scost2 = normalize_scost(scost2)
                
                mcosts2.append(mcost2)
                scosts2.append(scost2)

            pylab.figure()
            pylab.scatter(mcosts1, scosts1, c='b', label='greedy')
            pylab.plot(mcosts1, scosts1, c='b')

            pylab.scatter(mcosts2, scosts2, c='r', label='genetic')
            #pylab.plot(mcosts2, scosts2, c='r')

            pylab.legend()
            pylab.xlabel('spanning tree cost')
            pylab.ylabel('satellite cost')

            pylab.savefig('test_runtimes/pareto_front%d.pdf' % num_points)
            
            pylab.close()

            write_items = [num_points, greedy_runtime, genetic_runtime]
            
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
