import os

BUILDER_EXE = 'neuronbuilder2'

OUTDIR = 'neuron_builder'

def build_neuron(radius=1, output='output.swc'):
    build_command = '%s/%s %f > %s/%s' % (OUTDIR, BUILDER_EXE, radius,\
                                          OUTDIR, output)
    print build_command
    os.system(build_command)

def main():
    build_neuron()

if __name__ == '__main__':
    main()
