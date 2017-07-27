import networkx
import matplotlib as mpl
mpl.use('agg')
import pylab
from pareto_functions import pareto_prim, pareto_genetic
from time import time
from cost_functions import graph_costs
from neuron_utils import complete_graph, pareto_cost, sort_neighbors, is_tree
from random_graphs import random_point_graph
import argparse
import pandas as pd

MIN_POINTS = 8
MAX_POINTS = 50

COLNAMES = ['num_points', 'time1', 'time2', 'comparisons', 'dominates']

def runtimes_stats():
    df = pd.read_csv('test_runtimes.csv', names=COLNAMES)
    print "greedy success rate", float(sum(df['dominates'])) / float(sum(df['comparisons']))
    df = df[['num_points', 'time1', 'time2']]
    df = df.groupby('num_points', as_index=False).agg(pylab.mean)
    pylab.figure()
    time1 = df['time1'] * 60
    time2 = df['time2'] * 60
    pylab.plot(df['num_points'], time1, c='b', label='greedy')
    pylab.scatter(df['num_points'], time2, c='r', label='genetic')
    pylab.ylim(ymin=-1)
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
                mcost2, scost2, time2 = time_function(G, alpha, pareto_genetic)

                mcosts1.append(mcost1)
                mcosts2.append(mcost2)
                
                scosts1.append(scost1)
                scosts2.append(scost2)

                times1.append(time1)
                times2.append(time2)

                comparisons += 1
                cost1 = pareto_cost(mcost1, scost1, alpha)
                cost2 = pareto_cost(mcost2, scost2, alpha)
                if cost1 <= cost2:
                    dominates += 1

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

            time1 = pylab.mean(times1)
            time2 = pylab.mean(times2)

            write_items = [num_points, time1, time2, comparisons, dominates]
            
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
