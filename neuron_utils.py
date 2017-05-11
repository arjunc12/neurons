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
    print sorted(data)
    return data['species'], data['archive_name'], data['region1'], data['class1']

def main():
    neuron_file = argv[1]
    print map(str, neuron_info(neuron_file))

if __name__ == '__main__':
    main()
