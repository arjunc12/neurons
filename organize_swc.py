import os
from neuron_utils import neuron_info

SWC_DIR = 'neuron_nmo'
DATASET_DIR = 'datasets'

def remove_spaces(string):
    return string.replace(' ', '_')

def organize_swc():
    for filename in os.listdir(SWC_DIR):
        if filename[-4:] == '.swc':
            print filename
            info = neuron_info(filename)
            if info == None:
                print "nothing"
                continue
            print info
            cell_type, species, region, lab = map(remove_spaces, info)
            mv_dir = '%s/%s/%s/%s/%s' % (DATASET_DIR, cell_type, species, region, lab)
            os.system('mkdir -p %s' % mv_dir)
            os.system('mv %s/%s %s/%s' % (SWC_DIR, filename, mv_dir, filename))

def main():
    organize_swc()

if __name__ == '__main__':
    main()
