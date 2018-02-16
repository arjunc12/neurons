import argparse
import os
from collections import defaultdict

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'
PBS_DIR = '/iblsn/data/Arjun/neurons/pbs_files'
OUTPUT_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_output'
FIGS_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_figs'

def make_neuron_pbs(cell_types, animal_species, regions, labs, names, min_nodes, max_nodes):
    directory = DATASETS_DIR
    fpaths = defaultdict(list)
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
                        neuron_name = neuron_file[:-8]
                        if names != None and neuron_name not in names:
                            continue

                        if neuron_file[-8:] != ".CNG.swc": 
                            continue

                        bash_str = 'python pareto_steiner.py -n -c %s -s %s -r %s -l %s -na %s -min_nodes %d -max_nodes %d'\
                                                            % (cell_type, species,\
                                                               region, lab, neuron_name,\
                                                               min_nodes, max_nodes)
                        fname = 'pareto_steiner_%s.pbs' % (neuron_name)
                        fpath = '%s/%s' % (PBS_DIR, fname)
                        flines = []
                        if fpath not in fpaths:
                            flines += ['#!/bin/bash\n', 'cd /home/achandrasekhar/neurons/\n']
                        flines += '%s\n' % bash_str
                        fpaths[fpath] += flines
    
    for fpath in fpaths:
        flines = fpaths[fpath]
        with open(fpath, 'w') as f:
            f.writelines(flines)
        command = 'qsub %s' % fpath
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
