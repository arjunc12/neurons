import pandas as pd
from neuron_utils import DENDRITE_RATE

DATA_DRIVE = '/iblsn/data/Arjun/neurons'
CATEGORIES_FILE = '%s/neuron_categories/neuron_categories.csv' % DATA_DRIVE

OUTPUT_DIR = '%s/pareto_steiner_output' % DATA_DRIVE
PARETO_STEINER_FILE = '%s/pareto_steiner_synthetic%0.1f.csv' % (OUTPUT_DIR, DENDRITE_RATE)
MODELS_FILE = '%s/models_synthetic%0.1f.csv' % (OUTPUT_DIR, DENDRITE_RATE)
DENSITY_FILE = '%s/neuron_density/neuron_density.csv' % DATA_DRIVE
TRADEOFF_FILE = '%s/tradeoff_ratio_synthetic%0.1f.csv' % (OUTPUT_DIR, DENDRITE_RATE)
OUTPUT_FILE = '%s/neurons_df%0.1f.csv' % (OUTPUT_DIR, DENDRITE_RATE)

pareto_steiner_df = pd.read_csv(PARETO_STEINER_FILE, skipinitialspace=True)

categories_df = pd.read_csv(CATEGORIES_FILE, skipinitialspace=True)

models_df = pd.read_csv(MODELS_FILE, skipinitialspace=True)
models_df = models_df[models_df['model'] == 'neural']
models_df = models_df[['neuron name', 'neuron type', 'dist']]

density_df = pd.read_csv(DENSITY_FILE, skipinitialspace=True)
density_df['density'] = density_df['mcost'] = density_df['volume']
density_df = density_df[['neuron name', 'neuron type', 'mcost', 'volume', 'density']]

tradeoff_df = pd.read_csv(TRADEOFF_FILE, skipinitialspace=True)

df = pd.merge(pareto_steiner_df, categories_df)
df = pd.merge(df, models_df)
df = pd.merge(df, density_df)
df = pd.merge(df, tradeoff_df)

df.to_csv(OUTPUT_FILE, index=False)
