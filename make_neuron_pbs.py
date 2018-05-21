import argparse
import os
from collections import defaultdict
from random import sample

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'
PBS_DIR = '/iblsn/data/Arjun/neurons/pbs_files'
OUTPUT_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_output'
FIGS_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_figs'
DENSITIES_FILE = '/iblsn/data/Arjun/neurons/neuron_density/neuron_density.csv'
SIZE_FILE = '/iblsn/data/Arjun/neurons/neuron_size/neuron_size.csv'
SYNTHETIC_SIZE_FILE = '/iblsn/data/Arjun/neurons/neuron_size/neuron_size_synthetic.csv'
MAX_JOBS = 1000

NEURON_TYPES = {0 : 'axon', 1 : 'basal dendrite', 2 : 'apical dendrite',\
                3 : 'truncated axon'}

def read_sizes(size_file=SIZE_FILE):
    sizes = {}
    with open(size_file) as f:
        i = 0
        for line in f:
            i += 1
            if i == 1:
                continue
            line = line.strip('\n')
            line = line.split(', ')
            neuron_name = line[0]
            neuron_type = line[1]
            points = int(line[2])
            sizes[(neuron_name, neuron_type)] = points
    return sizes

def make_neuron_pbs(cell_types, animal_species, regions, labs, names,\
                    neuron_types, min_nodes, max_nodes, synthetic,\
                    max_jobs=MAX_JOBS):
    datasets_dir = DATASETS_DIR
    fpaths = defaultdict(list)
    size_file = None
    if synthetic:
        size_file = SYNTHETIC_SIZE_FILE
    else:
        size_file = SIZE_FILE
    neuron_sizes = read_sizes(size_file)
    for cell_type in os.listdir(datasets_dir):
        if cell_types != None and cell_type not in cell_types:
            continue
        for species in os.listdir(datasets_dir + '/' + cell_type):
            if animal_species != None and species not in animal_species:
                continue
            for region in os.listdir(datasets_dir + '/' + cell_type + '/' + species):
                if regions != None and region not in regions:
                    continue
                for lab in os.listdir(datasets_dir + "/" + cell_type + '/' + species+ '/' + region):
                    if labs != None and lab not in labs:
                        continue
                    for neuron_file in os.listdir(datasets_dir + "/" + cell_type + "/" + species + '/' + region + '/' + lab): 
                        neuron_name = neuron_file[:-8]
                        if names != None and neuron_name not in names:
                            continue

                        if neuron_file[-8:] != ".CNG.swc": 
                            continue

                        fname = 'pareto_steiner_%s.pbs' % (neuron_name)
                        fpath = '%s/%s' % (PBS_DIR, fname)
                        flines = []
                        
                    
                        for i in xrange(4):
                            if neuron_types != None and i not in neuron_types:
                                continue

                            neuron_type = NEURON_TYPES[i]
                            key = (neuron_name, neuron_type)
                            if key not in neuron_sizes:
                                continue
                            points = neuron_sizes[(neuron_name, neuron_type)]
                            if min_nodes <= points <= max_nodes:
                                if fpath not in fpaths:
                                    flines += ['#!/bin/bash\n', 'cd /home/achandrasekhar/neurons/\n']
                                bash_str = 'python pareto_steiner.py -n -c %s -s %s -r %s -l %s -na %s -min_nodes %d -max_nodes %d -nt %d'\
                                                                    % (cell_type,\
                                                                       species,\
                                                                       region, lab,\
                                                                       neuron_name,\
                                                                       min_nodes,\
                                                                       max_nodes,\
                                                                       i)
                                if synthetic:
                                    bash_str += ' --synthetic'
                                flines.append('%s\n' % bash_str)
                        
                        if len(flines) > 0:
                            fpaths[fpath] += flines

    print len(fpaths), "jobs"
    submit = raw_input('submit jobs? y/n: ')
    if submit != 'y':
        return None
    for fpath in fpaths:
        flines = fpaths[fpath]
        with open(fpath, 'w') as f:
            f.writelines(flines)
        command = 'qsub %s' % fpath
        print command
        os.system(command)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cell_types', nargs='+', default=None)
    parser.add_argument('-s', '-a', '--species', '--animal_species,', nargs='+',\
                        default=None, dest='animal_species')
    parser.add_argument('-r', '--regions', nargs='+', default=None)
    parser.add_argument('-l', '--labs', nargs='+', default=None)
    parser.add_argument('-na', '--names', nargs='+', default=None, dest='names')
    parser.add_argument('-nt', '--neuron_types', nargs='+', default=None)
    parser.add_argument('--synthetic', action='store_true')
    parser.add_argument('--min_nodes', type=int, default=50)
    parser.add_argument('--max_nodes', type=int, default=1000)
    parser.add_argument('--max_jobs', type=int, default=MAX_JOBS)

    args = parser.parse_args()

    cell_types = args.cell_types
    animal_species = args.animal_species
    regions = args.regions
    labs = args.labs
    names = args.names
    neuron_types = args.neuron_types
    if neuron_types != None:
        neuron_types = map(int, neuron_types)
    
    synthetic = args.synthetic
    
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes
    
    max_jobs = args.max_jobs

    make_neuron_pbs(cell_types, animal_species, regions, labs, names,\
                    neuron_types, min_nodes, max_nodes, synthetic, max_jobs)

if __name__ == '__main__':
    main()
