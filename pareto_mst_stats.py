from collections import defaultdict
import pandas as pd
from pareto_mst_plots import COLUMNS
import pylab

fname = 'pareto_mst.csv'
df = pd.read_csv(fname, names=COLUMNS)
total_trials = df['trials'].sum()
total_successes = df['successes'].sum()
print float(total_successes) / total_trials
print pylab.mean(df['neural_dist'] / df['centroid_dist'])
print pylab.mean(df['neural_dist'] / df['random_dist'])
