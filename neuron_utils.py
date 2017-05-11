import urllib2
import json
from sys import argv

def neuron_info(neuron_file):
    assert neuron_file[-4:] == '.swc'
    neuron = neuron_file[:-4]
    if neuron[-4:] == '.CNG':
        neuron = neuron[:-4]
    url = 'http://neuromorpho.org:8081/neuron/query/neuron_name&=&' + neuron
    data = urllib2.urlopen(url)
    data = json.load(data)
    data = data['data'][0]
    
    species = data['species']
    lab = data['archive_name']
    region = data['region1']
    class1 = data['class1']
    cell_type = class1
    if class1 == 'interneuron':
        class2 = data['class2']
        if class2 != None and class2 != '':
            cell_type = class2
    
    return cell_type, species, region, lab
    

def main():
    neuron_file = argv[1]
    print map(str, neuron_info(neuron_file))

if __name__ == '__main__':
    main()
