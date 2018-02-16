import os
import urllib2
import json
from sys import argv

SWC_DIR = 'neuron_nmo'
DATASET_DIR = '/iblsn/data/Arjun/neurons/datasets'

TEST_NEURON_NAME = 'p20'

def neuron_info(neuron_file):
    assert neuron_file[-4:] == '.swc'
    neuron = neuron_file[:-4]
    if neuron[-4:] == '.CNG':
        neuron = neuron[:-4]
    #url = 'http://neuromorpho.org:8081/neuron/query/neuron_name&=&' + neuron
    #url = 'http://neuromorpho.org/api/neuron/query/neuron_name/' + neuron
    url = 'http://neuromorpho.org/api/neuron/name/' + neuron
    data = urllib2.urlopen(url)
    data = json.load(data)
    if '_embedded' not in data:
        return None
    data = data['_embedded']
    if 'neuronResources' not in data:
        return None
    data = data['neuronResources']
    if len(data) == 0:
        return None
    data = data[0]
    
    if 'species' not in data:
        return None
    species = data['species']
    if species == None:
        return None
    if type(species) != list:
        species = [species]

    if 'archive' not in data:
        return None
    labs = data['archive']
    if labs == None:
        return None
    if type(labs) != list:
        labs = [labs]
    
    if 'brain_region' not in data:
        return None
    regions = data['brain_region']
    if regions == None:
        return None
    if type(regions) != list:
        regions = [regions]

    if 'cell_type' not in data:
        return None
    cell_types = data['cell_type']
    if cell_types == None:
        return None
    if type(cell_types) != list:
        cell_types = [cell_types]

    info_tuples = []
    for sps in species:
        if sps == None:
            continue
        for lab in labs:
            if lab == None:
                continue
            for region in regions:
                if region == None:
                    continue
                for cell_type in cell_types:
                    if cell_type == None:
                        continue
                    info_tuples.append((cell_type, sps, region, lab))

    return info_tuples

def remove_spaces(string):
    return string.replace(' ', '_')

def organize_swc():
    for filename in os.listdir(SWC_DIR):
        if filename[-4:] == '.swc':
            print filename
            info_tuples = neuron_info(filename)
            if info_tuples == None:
                print "nothing"
                continue
            for info in info_tuples:
                print info
                cell_type, species, region, lab = map(remove_spaces, info)
                mv_dir = '%s/%s/%s/%s/%s' % (DATASET_DIR, cell_type, species, region, lab)
                os.system('mkdir -p %s' % mv_dir)
                os.system('cp %s/%s %s/%s' % (SWC_DIR, filename, mv_dir, filename))

def test_api(neuron_name=TEST_NEURON_NAME):
    neuron_file = TEST_NEURON_NAME + '.swc'
    print neuron_info(neuron_file)

def main():
    if len(argv) > 1:
        for name in argv[1:]:
            test_api(name)
    else:
        organize_swc()

if __name__ == '__main__':
    main()
