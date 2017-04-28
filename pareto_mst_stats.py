from collections import defaultdict
import pandas as pd

columns = ['neuron', 'neural_dist', 'rand_dist']
fname = 'pareto_mst.csv'
df = pd.read_csv(fname, names=columns)

for name, group in df.groupby('neuron'):
    print name
    print group[group['rand_dist'] < group['neural_dist']].count()
