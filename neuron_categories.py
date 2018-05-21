import os
import pandas as pd

CATEGORIES_FILE = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories.csv'
CATEGORIES_FILE_FILTERED = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories_filtered.csv'
DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'

def remove_sep(string, sep='_'):
    return string.replace(sep, ' ')

def add_count_col(df, categories):
    return df.groupby(categories).size().reset_index(name='count')

def remove_commas(string):
    return string.replace(',', '')

def write_categories():
    with open(CATEGORIES_FILE, 'w') as f:
        f.write('neuron name, cell type, species, region, lab\n')
        for cell_type in os.listdir(DATASETS_DIR):
            for species in os.listdir(DATASETS_DIR + '/' + cell_type):
                for region in os.listdir(DATASETS_DIR + '/' + cell_type + '/' + species):
                    for lab in os.listdir(DATASETS_DIR + "/" + cell_type + '/' + species+ '/' + region):
                        for neuron_file in os.listdir(DATASETS_DIR + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                            if neuron_file[-8:] != ".CNG.swc": 
                                continue
                            name = neuron_file[:-8]
                            write_items = [name, cell_type.lower(), species.lower(), region.lower(), lab.lower()]
                            write_items = map(remove_commas, write_items)
                            write_items = ', '.join(write_items)
                            f.write('%s\n' % write_items)

def filter_categories():
    df = pd.read_csv(CATEGORIES_FILE, skipinitialspace=True)
    categories = ['cell type', 'species', 'region', 'lab']
    for category in ['cell type', 'species', 'region', 'lab']:
        df[category] = df[category].apply(remove_sep)
    df.drop_duplicates(inplace=True)
    
    for category in categories:
        df2 = add_count_col(df, category)
        df2 = pd.merge(df, df2)
        max_vals = {}
        for name, group in df2.groupby('neuron name'):
            counts = list(group['count'])
            vals = list(group[category])
            count_vals = zip(counts, vals)
            max_count, max_val = max(count_vals)
            max_vals[name] = max_val
        df['max ' + category] = df['neuron name'].map(max_vals)
        df = df[df[category] == df['max ' + category]]
        df.drop('max ' + category, axis=1, inplace=True)

    df.to_csv(CATEGORIES_FILE_FILTERED, index=False)

def main():
    #write_categories()
    filter_categories()

if __name__ == '__main__':
    main()
