import argparse
import os

def make_neuron_pbs(cell_types, animal_species, regions, labs, names, min_nodes, max_nodes):
    directory = 'datasets'
    fnames = []
    for cell_type in os.listdir(directory):
        if cell_types != None and cell_type not in cell_types:
            continue
        for species in os.listdir(directory + '/' + cell_type):
            if animal_species != None and species not in animal_species:
                continue
            for region in os.listdir(directory + '/' + cell_type + '/' + species):
                if regions != None and region not in regions:
                    continue
                for lab in os.listdir(directory + "/" + cell_type + '/' + species+ '/' + region):
                    if labs != None and lab not in labs:
                        continue
                    for neuron_file in os.listdir(directory + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = directory + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron_file
 
                        neuron_name = neuron_file[:-8]
                        if names != None and neuron_name not in names:
                            continue

                        if neuron_file[-8:] != ".CNG.swc": 
                            continue

                        bash_str = 'python pareto_steiner.py -n -c %s -s %s -r %s -l %s -na %s -min_nodes %d -max_nodes %d'\
                                                            % (cell_type, species,\
                                                               region, lab, neuron_name,\
                                                               min_nodes, max_nodes)
                        fname = 'pareto_steiner_%s.pbs' % neuron_name
                        fnames.append(fname)
                        f = open(fname, 'a')
                        f.write('#!/bin/bash\n')
                        f.write('cd $SGE_O_WORKDIR\n')
                        f.write('%s\n' % bash_str)
                        f.close()

    for fname in fnames:
        command = 'qsub %s' % fname
        print command
        os.system(command)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cell_types', nargs='+', default=None)
    parser.add_argument('-s', '-a', '--species', '--animal_species,', nargs='+',\
                        default=None, dest='animal_species')
    parser.add_argument('-na', '--names', nargs='+', default=None, dest='names')
    parser.add_argument('-r', '--regions', nargs='+', default=None)
    parser.add_argument('-l', '--labs', nargs='+', default=None)
    parser.add_argument('--min_nodes', type=int, default=50)
    parser.add_argument('--max_nodes', type=int, default=1000)

    args = parser.parse_args()

    cell_types = args.cell_types
    animal_species = args.animal_species
    regions = args.regions
    names = args.names
    labs = args.labs
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes

    make_neuron_pbs(cell_types, animal_species, regions, labs, names, min_nodes, max_nodes)

if __name__ == '__main__':
    main()
