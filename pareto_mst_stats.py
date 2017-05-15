from collections import defaultdict
import pandas as pd
from pareto_mst_plots import COLUMNS
import pylab
from numpy.ma import masked_invalid

def infmean(arr):
    return pylab.mean(masked_invalid(arr))

fname = 'pareto_mst.csv'
df = pd.read_csv(fname, names=COLUMNS)
total_trials = df['trials'].sum()
total_successes = df['successes'].sum()
print float(total_successes) / total_trials
print infmean(df['neural_dist'] / df['centroid_dist'])
print infmean(df['neural_dist'] / df['random_dist'])
