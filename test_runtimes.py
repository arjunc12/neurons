import networkx
import matplotlib as mpl
mpl.use('agg')
import pylab
from pareto_functions import pareto_prim, pareto_genetic
from time import time
from cost_functions import graph_costs
from neuron_utils import complete_graph, pareto_cost, sort_neighbors
from random_graphs import random_point_graph

MIN_POINTS = 8
MAX_POINTS = 15

def time_function(G, alpha, pareto_func):
    start = time()
    mst = pareto_func(G, alpha)
    end = time()
    runtime = (end - start) / 60.0
    mcost, scost = graph_costs(G)
    return mcost, scost, time

def test_runtimes(num_iters=10):
    for i in xrange(num_iters):
        print "iteration", i
        for num_points in xrange(MIN_POINTS, MAX_POINTS + 1):
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
            pylab.plot(mcosts2, scosts2, c='r')

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
    test_runtimes(num_iters=1)

if __name__ == '__main__':
    main()
