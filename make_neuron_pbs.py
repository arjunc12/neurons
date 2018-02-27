import argparse
import os
from collections import defaultdict
from random import sample

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'
PBS_DIR = '/iblsn/data/Arjun/neurons/pbs_files'
OUTPUT_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_output'
FIGS_DIR = '/iblsn/data/Arjun/neurons/pareto_steiner_output/steiner_figs'
MAX_JOBS = 1000

def make_neuron_pbs(cell_types, animal_species, regions, labs, names,\
                    min_nodes, max_nodes, max_jobs=MAX_JOBS):
    datasets_dir = DATASETS_DIR
    fpaths = defaultdict(list)
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
                        
                        if fpath not in fpaths:
                            flines += ['#!/bin/bash\n', 'cd /home/achandrasekhar/neurons/\n']
                    
                        bash_str = 'python pareto_steiner.py -n -c %s\
                                                             -s %s -r %s\
                                                             -l %s -na %s\
                                                             -min_nodes %d\
                                                             -max_nodes %d'\
                                                            % (cell_type,\
                                                               species,\
                                                               region, lab,\
                                                               neuron_name,\
                                                               min_nodes,\
                                                               max_nodes)
                        flines += '%s\n' % bash_str
                        
                        fpaths[fpath] += flines

    print len(fpaths)
    sample_size = min(len(fpaths.keys()), max_jobs)
    for fpath in sample(fpaths.keys(), sample_size):
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
    parser.add_argument('--max_jobs', type=int, default=MAX_JOBS)

    args = parser.parse_args()

    cell_types = args.cell_types
    animal_species = args.animal_species
    regions = args.regions
    names = args.names
    labs = args.labs
    min_nodes = args.min_nodes
    max_nodes = args.max_nodes
    max_jobs = args.max_jobs

    make_neuron_pbs(cell_types, animal_species, regions, labs, names,\
                    min_nodes, max_nodes, max_jobs)

if __name__ == '__main__':
    main()
