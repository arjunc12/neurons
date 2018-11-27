from pareto_steiner import *
from neuron_utils import *
import os
from collections import defaultdict
from random import shuffle

RATE = 4

MIN_COUNT = 100
MIN_POINTS = 250
MAX_POINTS = 750

COMPARISON_FRONTS_DIR = '/iblsn/data/Arjun/neurons/neuron_types_comparison%0.1f/fronts' % RATE
COMPARISON_OUTPUT_DIR = '/iblsn/data/Arjun/neurons/neuron_types_comparison%0.1f/output' % RATE
COMPARISON_FIGS_DIR = '/iblsn/data/Arjun/neurons/neuron_types_comparison%0.1f/figs' % RATE

output = True

datasets_dir = NEUROMORPHO_DATASETS_DIR


def main():
    neuron_type_counts = {}
    neuron_type_counts['axon'] = 0
    neuron_type_counts['apical dendrite'] = 0
    neuron_type_counts['basal dendrite'] = 0
    cell_types = os.listdir(datasets_dir)
    shuffle(cell_types)
    for cell_type in cell_types:
        if cell_type in SKIP_TYPES:
            continue
        specie = os.listdir(datasets_dir + '/' + cell_type)
        shuffle(specie)
        for species in specie:
            regions = os.listdir(datasets_dir + '/' + cell_type + '/' + species)
            shuffle(regions)
            for region in regions:
                labs = os.listdir(datasets_dir + "/" + cell_type + '/' + species+ '/' + region)
                shuffle(labs)
                for lab in labs:
                    for neuron_file in os.listdir(datasets_dir + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = datasets_dir + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron_file
                        
                        neuron_name = neuron_file[:-8]

                        if neuron_file[-8:] != ".CNG.swc": 
                            continue
                        
                        try:
                            graphs = get_neuron_points(filename)
                        except AssertionError as e:
                            continue

                        for i, G in enumerate(graphs):
                            if G == None:
                                continue
 
                            neuron_type = NEURON_TYPES[i]
                            if neuron_type == 'truncated axon':
                                continue
                            if neuron_type_counts[neuron_type] >= MIN_COUNT:
                                continue
                            ntype = neuron_type.replace(' ', '_')
                            
                            H = add_synapses(G, neuron_type=neuron_type, rate=RATE)
                            if not MIN_POINTS <= len(H.graph['synapses']) + 1 <= MAX_POINTS:
                                continue

                            fronts_dir = COMPARISON_FRONTS_DIR + '/' + neuron_name + '/' + ntype
                            fronts_dir = fronts_dir.replace(' ', '_')
                            output_dir = COMPARISON_OUTPUT_DIR
                            figs_dir = COMPARISON_FIGS_DIR + '/' + neuron_name + '/' + ntype
                            
                            for directory in [fronts_dir, output_dir]:
                                os.system('mkdir -p %s' % directory)

                            pareto_analysis(H, neuron_name, neuron_type,\
                                            fronts_dir=fronts_dir,\
                                            output_dir=output_dir,\
                                            figs_dir=figs_dir, output=True,\
                                            viz_trees=False)

                            neuron_type_counts[neuron_type] += 1

                            if min(neuron_type_counts.values()) == MIN_COUNT:
                                return None

if __name__ == '__main__':
    main()
