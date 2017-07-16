import os
from neuron_utils import get_neuron_points, viz_tree

BUILDER_EXE = 'neuronbuilder2'

OUTDIR = 'neuron_builder'
OUTFILE = 'neuron_builder.swc'

def build_neuron(radius=1):
    build_command = './%s %f > %s' % (BUILDER_EXE, radius, OUTFILE)
    print build_command
    os.system(build_command)
    G = get_neuron_points(OUTFILE)[0]
    return G

def main():
    G = build_neuron(1)
    viz_tree(G, 'neuron_builder', OUTDIR)

if __name__ == '__main__':
    main()
