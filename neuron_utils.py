import matplotlib as mpl
mpl.use('agg')
import urllib2
import json
from sys import argv
import networkx as nx
import pylab
import numpy as np

def centroid(G):
    root = G.graph['root']
    root_coord = G.node[root]['coord']
    centroid = np.zeros(len(root_coord))
    for u in G.nodes_iter():
        point = G.node[u]['coord']
        assert len(point) == len(root_coord)
        if point != root:
            centroid += point
    centroid /= G.number_of_nodes() - 1
    return centroid

def label_points(G):
    root = G.graph['root']
    for u in G:
        if u == root: 
            G.node[u]['label'] = "root"
        elif G.degree(u) == 1: 
            G.node[u]['label'] = "tip"
        elif G.degree(u) > 2: 
            G.node[u]['label'] = "branch"
        else:                  
            G.node[u]['label'] = "continue"

def viz_tree(G, name, outdir='figs'):
    """ Displays plant/tree visualization. """
    
    root = G.graph['root']
    node_size,node_color = [],[]
    pos = {}    
    for u in G:
        coord = G.node[u]['coord']
        pos[u] = (coord[0], coord[1])
        
        label = G.node[u]['label']
        if label == "root":
            node_color.append('black')
            node_size.append(350)
        elif label in ['branch', 'tip']:
            node_color.append('green')
            node_size.append(100)
        elif label == "continue":
            node_color.append('brown')
            #node_size.append(250)   
            node_size.append(1)   
        elif label == 'centroid':
            node_color.append('purple')
            node_size.append(250)
        elif label == 'isolated_start':
            node_color.append('yellow')
            node_size.append(100)
        elif label == 'isolated_end':
            node_color.append('blue')
            node_size.append(100)
        else:
            assert False

    nx.draw(G,pos=pos,arrows=False,with_labels=False,node_size=node_size,node_color=node_color,edge_color="brown",width=4,font_size=12,font_color='red',font_weight='bold')
    #nx.draw(G,pos=pos,arrows=False,with_labels=False,node_size=15,font_size=12,font_color='red',font_weight='bold')
    pylab.draw()
    #PP.show()
    
    '''
    mcost = mst_cost(G)
    scost = satellite_cost(G)
    title_str = 'satellite cost = %f\nspanning tree cost = %f' % (scost, mcost)
    #pylab.text(125, 125, title_str)
    '''

    pylab.savefig("%s/%s.pdf" % (outdir, name))
    pylab.close()

def point_dist(p1, p2):
    assert len(p1) == len(p2)
    sq_dist = 0
    for i in xrange(len(p1)):
        x1, x2 = p1[i], p2[i]
        sq_dist += (x1 - x2) ** 2
    return sq_dist ** 0.5

def neuron_info(neuron_file):
    assert neuron_file[-4:] == '.swc'
    neuron = neuron_file[:-4]
    if neuron[-4:] == '.CNG':
        neuron = neuron[:-4]
    url = 'http://neuromorpho.org:8081/neuron/query/neuron_name&=&' + neuron
    data = urllib2.urlopen(url)
    data = json.load(data)
    data = data['data']
    if len(data) == 0:
        return None
    data = data[0]
    
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
