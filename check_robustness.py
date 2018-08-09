import pandas as pd

CELL_TYPES = ['interneuron', 'sensory receptor', 'pyramidal', 'basket', 'martinotti',\
              'neurogliaform', 'fast-spiking', 'golgi', 'chandelier',\
              'stellate', 'thick-tufted', 'purkinje', 'granule', 'mitral',\
              'amacrine', 'projection', 'serotonergic', 'gabaergic',\
              'nitrergic', 'cholinergic', 'glutamatergic']

df = pd.read_csv('neurons_df0.2.csv', skipinitialspace=True)
df = df[df['points'] >= 0]

for cell_type in CELL_TYPES:
    print "---------------------"
    print cell_type
    df2 = df[df['cell type'] == cell_type]
    print len(set(zip(df2['neuron name'], df2['neuron type']))), 'unique neurons'
    print len(df2['species'].unique()), "unique species"
