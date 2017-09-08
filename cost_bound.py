from random_graphs import random_point_graph
from pareto_functions import pareto_prim
import pylab
from cost_functions import best_mst_cost, best_satellite_cost, graph_costs
import argparse
import pandas as pd

BOUND = 100

def cost_bounds(points):
    print "points", points
    G = random_point_graph(points, -BOUND, BOUND, -BOUND, BOUND, -BOUND, BOUND)
    delta = 0.01
    alphas = pylab.arange(delta, 1, delta)
    opt_mcost = best_mst_cost(G)
    opt_scost = best_satellite_cost(G)
    f = open('cost_bound.csv', 'a')
    for alpha in alphas:
        print alpha
        pareto_mst = pareto_prim(G, alpha)
        mcost, scost = graph_costs(pareto_mst)
        
        write_items = [points, alpha, mcost, scost, opt_mcost, opt_scost]
        write_items = map(str, write_items)
        write_items = ', '.join(write_items)
        f.write('%s\n' % write_items)

def cost_bound_stats():
    df = pd.read_csv('cost_bound.csv', names=['points', 'alpha', 'mcost',\
                                              'scost', 'opt_mcost', 'opt_scost'],\
                                              skipinitialspace=True)
    df['mcost_ratio'] = df['mcost'] / df['opt_mcost']
    df['scost_ratio'] = df['scost'] / df['opt_scost']

    print "max mcost ratio", max(df['mcost_ratio'])
    print "max scost ratio", max(df['scost_ratio'])

    max_mcost_ratios = []
    max_scost_ratios = []
    points = []
    for name, group in df.groupby('points'):
        points.append(name)
        max_mcost_ratios.append(max(group['mcost_ratio']))
        max_scost_ratios.append(max(group['scost_ratio']))
 
    pylab.figure()
    pylab.scatter(points, max_mcost_ratios)
    pylab.savefig('cost_bounds/mcost_bounds.pdf', format='pdf')
    pylab.close()

    pylab.figure()
    pylab.scatter(points, max_scost_ratios)
    pylab.savefig('cost_bounds/scost_bounds.pdf', format='pdf')
    pylab.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--min_points', type=int, default=100)
    parser.add_argument('--max_points', type=int, default=1000)
    parser.add_argument('-x', '--num_iters', type=int, default=10)
    parser.add_argument('-b', '--cost_bounds', action='store_true')
    parser.add_argument('-s', '--stats', action='store_true')

    args = parser.parse_args()
    min_points = args.min_points
    max_points = args.max_points
    num_iters = args.num_iters
    bounds = args.cost_bounds
    stats = args.stats

    if stats:
        cost_bound_stats()

    if bounds:
        for i in xrange(num_iters):
            for points in xrange(min_points, max_points + 1):
                cost_bounds(points)

if __name__ == '__main__':
    main()
