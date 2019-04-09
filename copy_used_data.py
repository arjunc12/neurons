import pandas as pd
from pareto_steiner_stats import get_dfs
import os

DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'

output = pd.read_csv('/iblsn/data/Arjun/neurons/pareto_steiner_output0.2/pareto_steiner_synthetic0.2.csv', skipinitialspace=True)
output = output[output['neuron type'] != 'truncated axon']
output = output[output['points'] >= 100]

categories = pd.read_csv('/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories.csv', skipinitialspace=True)

df = pd.merge(output, categories)
unique_names = set(df['neuron name'].unique())

for cell_type in os.listdir(DATASETS_DIR):
    for species in os.listdir('/'.join([DATASETS_DIR, cell_type])):
        for region in os.listdir('/'.join([DATASETS_DIR, cell_type, species])):
            for lab in os.listdir('/'.join([DATASETS_DIR, cell_type, species, region])):
                for fname in os.listdir('/'.join([DATASETS_DIR, cell_type, species, region, lab])):
                    name = fname[:-8]
                    if name in unique_names:
                        print name
                        os.system('cp %s/%s/%s/%s/%s/%s tracings_used/' % (DATASETS_DIR, cell_type, species, region, lab, fname))
