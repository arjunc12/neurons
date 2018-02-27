import os

CATEGORIES_FILE = '/iblsn/data/Arjun/neurons/neuron_categories/neuron_categories.csv'
DATASETS_DIR = '/iblsn/data/Arjun/neurons/datasets'

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

def main():
    write_categories()

if __name__ == '__main__':
    main()
