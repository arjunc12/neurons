import pandas as pd
NEURONS_FILE = '/iblsn/data/Arjun/neurons/pareto_steiner_output/neurons_df0.2.csv'

INTERESTING_CELL_TYPES = ['interneuron', 'sensory receptor', 'pyramidal',\
                          'basket', 'martinotti', 'neurogliaform',\
                          'fast-spiking', 'golgi', 'chandelier', 'stellate',\
                          'purkinje', 'granule', 'mitral', 'amacrine',\
                          'projection'] 

INTERESTING_TRANSMITTERS = ['serotonergic', 'gabaergic', 'nitrergic',\
                            'cholinergic', 'glutamatergic']

CELL_TYPES = INTERESTING_CELL_TYPES + INTERESTING_TRANSMITTERS

def main():
    df = pd.read_csv(NEURONS_FILE, skipinitialspace=True)
    df = df[df['points'] >= 0]

    for cell_type in CELL_TYPES:
        print "---------------------"
        print cell_type
        df2 = df[df['cell type'] == cell_type]
        print len(set(zip(df2['neuron name'], df2['neuron type']))), 'unique neurons'
        print len(df2['species'].unique()), "unique species"

if __name__ == '__main__':
    main()
