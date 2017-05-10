import networkx as nx
import matplotlib as mpl
mpl.use('agg')
import pylab
from pareto_mst import *
from sys import argv

def draw_spanning_tree(G, P2Coord, P2Label, name):
    mst = min_spanning_tree(G)
    viz_tree(mst, P2Coord, P2Label, name + '_mst')

def draw_satellite_tree(G, root, P2Coord, P2Label, name):
    satellite = satellite_tree(G, root)
    viz_tree(satellite, P2Coord, P2Label, name + '_satellite')

def draw_neuron(neuron_file):
    G, P2Coord, root = get_neuron_points(neuron_file)
    P2Label = label_points(G, root)
    assert neuron_file[-4:] == '.swc'
    name = neuron_file[9:-4]
    viz_tree(G, P2Coord, P2Label, name)
    initialize_lengths(G)
    draw_spanning_tree(G, P2Coord, P2Label, name)
    draw_satellite_tree(G, root, P2Coord, P2Label, name)

def draw_pareto_trees(neuron_file):
    delta = 0.01
    alphas = pylab.arange(0, delta, 1 + delta)
    for alpha in alphas:
        pass

def main():
    neuron_file = argv[1]
    draw_neuron(neuron_file)

if __name__ == '__main__':
    main()
