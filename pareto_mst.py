import numpy as np
import networkx as nx
from sys import argv
import matplotlib as mpl
mpl.use('agg')
import pylab
import os
from random import shuffle
from itertools import combinations

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
        else:
            assert False

    nx.draw(G,pos=pos,arrows=False,with_labels=False,node_size=node_size,node_color=node_color,edge_color="brown",width=4,font_size=12,font_color='red',font_weight='bold')
    #nx.draw(G,pos=pos,arrows=False,with_labels=False,node_size=15,font_size=12,font_color='red',font_weight='bold')
    pylab.draw()
    #PP.show()

    mcost = mst_cost(G)
    scost = satellite_cost(G, root)
    title_str = 'satellite cost = %f\nspanning tree cost = %f' % (scost, mcost)
    pylab.text(125, 125, title_str)

    pylab.savefig("%s/%s.pdf" % (outdir, name))
    pylab.close()

def get_neuron_points(filename, dim='3D'):
    #for arbor_type in ["2","3","4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
    H = nx.Graph()
    H_root = None
    HP2Coord = None
    for arbor_type in ["2", "3", "4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
        
        # Reads in 3D arborization.
        G = nx.Graph()
        P2Coord = {} # u-> (x,y)
        root = -1
        with open(filename) as f:
            for line in f:
                if line.startswith("#"): continue

                cols = line.strip().split()
                assert len(cols) == 7

                if not (cols[1] == "1" or cols[1] == arbor_type): continue
                
                if cols[6] == "-1":
                    if root != -1: assert False # duplicate root.

                    root = int(cols[0])
                                        
                    assert root not in P2Coord and root not in G
                    G.add_node(root)
                    coord = None
                    if dim == "3D":
                        coord = (float(cols[2]),float(cols[3]),float(cols[4]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        coord = (float(cols[2]),float(cols[3]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]))
                    else:
                        assert False
                    G.node[root]['coord'] = coord
                    G.graph['root'] = root


                else:
                    u,v = int(cols[6]),int(cols[0])
                    assert u in G and u in P2Coord
                    assert not G.has_edge(u,v)
            
                    coord = None
                    if dim == "3D":
                        coord = (float(cols[2]),float(cols[3]),float(cols[4]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        coord = (float(cols[2]),float(cols[3]))
                        #P2Coord[root] = (float(cols[2]),float(cols[3]))
                    else:
                        assert False
                    
                    G.node[v]['coord'] = coord
                    
                    G.add_edge(u, v)
                    G[u][v]['length'] = point_dist(G.node[u]['coord'], G.node[v]['coord'])

        if G.number_of_nodes() > H.number_of_nodes():
            H = G

    label_points(H)

    return H

def initialize_lengths(G):
    for u, v in G.edges_iter():
        G[u][v]['length'] = point_dist(u, v)

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

def makes_cycle(u, v, node_to_forest):
    f1 = node_to_forest[u]
    f2 = node_to_forest[v]
    return f1 == f2
 
def combine_forests(u, v, node_to_forest, forest_to_nodes, forest_to_edges):
    f1 = node_to_forest[u]
    f2 = node_to_forest[v]
    if f1 == f2:
        return
    new_forest = min(f1, f2)
    old_forest = max(f1, f2)
    update_nodes = forest_to_nodes[old_forest]
    update_edges = forest_to_edges[old_forest]
    for node in update_nodes:
        node_to_forest[node] = new_forest
    forest_to_nodes[new_forest] += update_nodes
    forest_to_edges[new_forest] += update_edges
    forest_to_edges[new_forest].append((u, v))
    del forest_to_nodes[old_forest]
    del forest_to_edges[old_forest]
    
def kruskal(nodes, edges):
    node_to_forest = {}
    forest_to_nodes = {}
    forest_to_edges = {}
    
    for i, u in enumerate(nodes):
        node_to_forest[u] = i
        forest_to_nodes[i] = [u]
        forest_to_edges[i] = []
        
    for u, v in edges:
        if not makes_cycle(u, v, node_to_forest):
            combine_forests(u, v, node_to_forest, forest_to_nodes, forest_to_edges)
            
    return sum(forest_to_edges.values(), [])

def random_mst(G):
    rand_edges = G.edges()
    shuffle(rand_edges)
    mst_edges = kruskal(G.nodes(), rand_edges)
    mst = nx.Graph()
    for u, v in mst_edges:
        mst.add_edge(u, v)
        mst[u][v]['length'] = point_dist(u, v)
    return mst

def mst_cost(G):
    total_length = 0
    for u, v in G.edges_iter():
        total_length += G[u][v]['length']
    return total_length

def satellite_cost(G):
    root = G.graph['root']
    total_cost = 0
    for u in G.nodes_iter():
        if u == root:
            continue
        if P2Label != None and P2Label[u] == 'continue':
            continue
        if nx.has_path(G, root, u):
            total_cost += nx.shortest_path_length(G, root, u, weight='length')
        else:
            total_cost += float("inf")
    return total_cost

def pareto_kruskal(G, alpha):
    root = G.graph['root']
    pareto_mst = nx.Graph()
    H = G.copy()
   
    pareto_mst.add_node(root)
    pareto_mst.node[root]['droot'] = 0
   
    graph_mcost = 0
    graph_scost = 0

    print "sorting neighbors"
    closest_neighbors = {}
    for u in G.nodes_iter():
        closest_neighbors[u] = sorted(G.neighbors(u), key=lambda v : G[u][v]['length'])

    while pareto_mst.number_of_nodes() < H.number_of_nodes():
        best_edge = None
        best_cost = float("inf")

        candidate_edges = set()
        for u in pareto_mst.nodes():
            assert 'droot' in pareto_mst.node[u]
             
            invalid_neighbors = []
            closest_neighbor = None
            for i in xrange(len(closest_neighbors[u])):
                v = closest_neighbors[u][i]
                if pareto_mst.has_node(v):
                    invalid_neighbors.append(v)
                else:
                    closest_neighbor = v
                    break

            for n in invalid_neighbors:
                closest_neighbors[u].remove(n)

            if closest_neighbor != None:
                candidate_edges.add(tuple(sorted((u, closest_neighbor))))

        for i, (u, v) in enumerate(candidate_edges):

            scost = graph_scost           
            if pareto_mst.has_node(u):
                scost += pareto_mst.node[u]['droot'] + G[u][v]['length']
            elif pareto_mst.has_node(v):
                scost += pareto_mst.node[v]['droot'] + G[u][v]['length']
            else:
                raise ValueError('something is wrong')

            mcost = graph_mcost + G[u][v]['length']
 
            mcost *= alpha
            scost *= (1 - alpha)
            cost = mcost + scost
            
            if cost < best_cost:
                best_edge = (u, v)
                best_cost = cost
        if best_edge == None:
            break
        u, v = best_edge
        pareto_mst.add_edge(u, v)
        pareto_mst[u][v]['length'] = G[u][v]['length']
        

        if 'droot' in pareto_mst.node[u]:
            pareto_mst.node[v]['droot'] = pareto_mst.node[u]['droot'] + pareto_mst[u][v]['length']
        elif 'droot' in pareto_mst.node[v]:
            pareto_mst.node[u]['droot'] = pareto_mst.node[v]['droot'] + pareto_mst[u][v]['length']
        else:
            raise ValueError('something is really wrong')

        H.remove_edge(u, v)
        closest_neighbors[u].remove(v)
        closest_neighbors[v].append(u)
    
    return pareto_mst 

def point_dist(p1, p2):
    assert len(p1) == len(p2)
    sq_dist = 0
    for i in xrange(len(p1)):
        x1, x2 = p1[i], p2[i]
        sq_dist += (x1 - x2) ** 2
    return sq_dist ** 0.5

def centroid(G):
    root = G.graph['root']
    centroid = np.zeros(len(root))
    for u in G.nodes_iter():
        point = G.node[u]['coord']
        assert len(point) == len(root)
        if point != root:
            centroid += point
    centroid /= len(points) - 1
    return centroid

def centroid_mst(G):
    cent_mst = G.copy()
    cent_mst.remove_edges_from(G.edges())
    
    centroidp = centroid(G)
    cen_mst.add_node(centroidp)
    cent_mst.node[centroidp]['coord'] = centroidp
    for u in G.nodes_iter():
        cent_mst.add_edge(u, centroidp)
        cent_mst[u][centroidp]['length'] = point_dist(cent_mst.node[u]['coord'], centroidp)
    return cent_mst

def centroid_mst_costs(G):
    root = G.graph['root']
    centroidp = centroid(G)
    root_cdist = point_dist(root, centroidp)

    mcost = root_cdist
    scost = 0
    
    for point in points:
        if point != root:
            cdist = point_dist(point, centroidp)
            mcost += cdist
            scost += cdist + root_cdist

    return mcost, scost

def points_to_graph(points):
    point_graph = nx.Graph()
    for p1, p2 in combinations(points, 2)
        point_graph.add_edge(p1, p2)
        length = point_dist(p1, p2)
        point_graph[p1][p2]['length'] = length
    return point_graph

def complete_graph(G):
    H = nx.copy(G)
    H.remove_edges_from(H.edges())
    for u, v in combinations(H.nodes(), 2):
        H.add_edge(u, v)
        H[u][v] = point_dist(H.node[u]['coord'], H.node[v]['coord'])
    return H

def pareto_dist(pareto_mcosts, pareto_scosts, mcost, scost):
    best_dist = float("inf")
    best_index = None

    assert len(pareto_mcosts) == len(pareto_scosts)

    for i in xrange(len(pareto_mcosts)):
        pareto_mcost = pareto_mcosts[i]
        pareto_scost = pareto_scosts[i]

        dist = point_dist((pareto_mcost, pareto_scost), (mcost, scost))
        if dist < best_dist:
            best_dist = dist
            best_index = i

    return best_dist, best_index

def satellite_tree(G):
    root = G.graph['root']
    
    satellite = G.copy()
    satellite.remove_edges_from(G.edges())

    for u in satellite.nodes_iter():
        if u != root:
            satellite.add_edge(u, root)
            p1, p2 = satellite.node[u]['coord'], satellite.node[root]['coord']
            satellite[u][root]['length'] = point_dist(p1, p2)
    
    return satellite

def best_satellite_cost(G):
    best_cost = 0
    root = G.graph['root']
    root_coord = G.node[root]['coord']
    for u in G.nodes_iter():
        best_cost += point_dist(root_coord, G.node[u]['coord'])
    return best_cost

def min_spanning_tree(G):
    sorted_edges = sorted(G.edges(), key= lambda x : G[x[0]][x[1]]['length'])
    mst_edges = kruskal(G.nodes(), sorted_edges)
    mst = G.copy()
    mst.remove_edges_from(mst.edges())
    mst.add_edges_from(mst_edges)
    return mst

def best_mst_cost(G):
    sorted_edges = sorted(G.edges(), key= lambda x : G[x[0]][x[1]]['length'])
    mst_edges = kruskal(G.nodes(), sorted_edges)
    best_cost = 0
    for u, v in mst_edges:
        best_cost += G[u][v]['length']
    return best_cost 

def non_continue_nodes(G, root):
    points = []
    for u in G.nodes_iter():
        degree = len(G.neighbors(u))
        if (u == root) or (u != root and degree != 2):
            points.append(u)
    return points

def non_continue_subgraph(G, root, P2Coord):
    non_cont_nodes = non_continue_nodes(G, root)
    H = nx.Graph()
    for u, v in combinations(non_cont_nodes, 2):
        H.add_edge(u, v)
        H[u][v]['length'] = point_dist(P2Coord[u], P2Coord[v])
    return H

def graph_to_points(G, root, P2Coord, cont=False):
    points = []
    for u in G.nodes_iter():
        degree = len(G.neighbors(u))
        if (u == root) or (u != root and degree != 2) or cont:
            points.append(P2Coord[u])
    return points

def init_lengths(G, P2Coord):
    for u, v in G.edges_iter():
        p1, p2 = P2Coord[u], P2Coord[v]
        G[u][v]['length'] = point_dist(p1, p2)

def pareto_drawings(filename, name, outdir='drawings'): 
    G, P2Coord, root_node = get_neuron_points(filename)
    P2Label = label_points(G, root_node)
    H = non_continue_subgraph(G, root_node, P2Coord)
    if H.number_of_nodes() <= 0 or H.number_of_nodes() > 1000:
        return None
 
    outdir = '%s/%s' % (outdir, name)
    os.system('mkdir -p %s' % outdir)

    delta = 0.01
    alphas = np.arange(0, 1 + delta, delta)
    for alpha in alphas:
        print "alpha", alpha
        pareto_tree = pareto_kruskal(H, root_node, alpha)
        viz_tree(pareto_tree, P2Coord, P2Label, name + str(alpha), outdir=outdir)


    init_lengths(G, P2Coord)
    viz_tree(G, P2Coord, P2Label, name + '_neural', outdir=outdir)

    sat_tree = satellite_tree(H, root_node)
    init_lengths(sat_tree, P2Coord)
    viz_tree(sat_tree, P2Coord, P2Label, name + '_satellite', outdir=outdir)

    mst = min_spanning_tree(H)
    init_lengths(mst, P2Coord)
    viz_tree(mst, P2Coord, P2Label, name + '_mst', outdir=outdir)

def pareto_plot(filename, name, cell_type, species, region, lab, outdir=None):
    G, P2Coord, root_node = get_neuron_points(filename)
    #points = P2Coord.values()
    points = graph_to_points(G, root_node, P2Coord)
    root = P2Coord[root_node]
    print len(points), "points"
    print G.number_of_nodes(), "nodes"
    if G.number_of_nodes() <= 0 or G.number_of_nodes() > 1000:
        return None
   
    print "making graph"
    point_graph = points_to_graph(points)
    
    mcosts = []
    scosts = []
    delta = 0.01
    alphas = np.arange(0, 1 + delta, delta)
    for alpha in alphas:
        print "alpha", alpha
        mst = pareto_kruskal(point_graph, root, alpha)
        mcost = mst_cost(mst)
        scost = satellite_cost(mst, root)
        
        mcosts.append(mcost)
        scosts.append(scost)

    pylab.plot(mcosts, scosts, c = 'b')
    pylab.scatter(mcosts, scosts, c='b', label='pareto mst')
    pylab.xlabel('spanning tree cost')
    pylab.ylabel('satellite cost')
    neural_mst = nx.Graph()
    for u, v in G.edges_iter():
        p1, p2 = P2Coord[u], P2Coord[v]
        neural_mst.add_edge(p1, p2)
        neural_mst[p1][p2]['length'] = point_dist(p1, p2)

    neural_mcost = mst_cost(neural_mst)
    neural_scost  = satellite_cost(neural_mst, root)

    neural_dist, neural_index = pareto_dist(mcosts, scosts, neural_mcost, neural_scost)
    neural_closem = mcosts[neural_index] 
    neural_closes = scosts[neural_index]
    neural_alpha = alphas[neural_index]

    pylab.scatter([neural_mcost], [neural_scost], c='r', marker='x', linewidths=15, label='neural mst')
    pylab.plot([neural_mcost, neural_closem], [neural_scost, neural_closes], c='r', linestyle='--')

    centroid_mcost, centroid_scost = centroid_mst_costs(points, root)
    centroid_dist, centroid_index = pareto_dist(mcosts, scosts, centroid_mcost, centroid_scost)
    centroid_closem = mcosts[centroid_index]
    centroid_closes = scosts[centroid_index]
    centroid_alpha = alphas[centroid_index]

    pylab.scatter([centroid_mcost], [centroid_scost], c='g', marker='+', linewidths=15, label='centroid mst')
    pylab.plot([centroid_mcost, centroid_closem], [centroid_scost, centroid_closes], c='g', linestyle='--')

    opt_scost = best_satellite_cost(root, points)
    opt_mcost = best_mst_cost(point_graph)

    #pylab.axhline(opt_scost)
    #pylab.axvline(opt_mcost)

    if outdir == None:
        outdir = 'figs'

    pylab.legend()

    pylab.savefig('%s/pareto_mst_%s.pdf' % (outdir, name), format='pdf')
    pylab.close()

    f = open('pareto_mst.csv', 'a')
    ntrials = 100
    successes = 0
    total_rand_dist = 0.0
    for i in xrange(ntrials):
        rand_mst = random_mst(point_graph)
        rand_mcost = mst_cost(rand_mst)
        rand_scost = satellite_cost(rand_mst, root)
        rand_dist, rand_index = pareto_dist(mcosts, scosts, rand_mcost, rand_scost)
        total_rand_dist += rand_dist
        rand_closem = mcosts[rand_index]
        rand_closes = scosts[rand_index]
        rand_alpha = alphas[rand_index]
        if rand_dist < neural_dist:
            successes += 1
    mean_rand_dist = total_rand_dist / ntrials
    write_items = [name, cell_type, species, region, lab]
    write_items.append(neural_alpha)
    write_items += [neural_dist, centroid_dist, mean_rand_dist]
    write_items += [ntrials, successes]
    write_items = map(str, write_items)
    write_items = ', '.join(write_items)
    f.write('%s\n' % write_items)

    #pylab.scatter([rand_mcost], [rand_scost], c='m', marker='*', linewidths=15)
    #pylab.plot([rand_mcost, rand_closem], [rand_scost, rand_closes], c='m', linestyle='-')


def neuromorpho_plots(plot_species=None):
    #directory = 'neuromorpho'
    directory = 'datasets'
    i = 0
    for cell_type in os.listdir(directory):
        for species in os.listdir(directory + '/' + cell_type):
            if plot_species != None and species not in plot_species:
                continue
            for region in os.listdir(directory + '/' + cell_type + '/' + species):
                for lab in os.listdir(directory + "/" + cell_type + '/' + species+ '/' + region):
                    for neuron in os.listdir(directory + "/" + cell_type + "/" + species + '/' + region + '/' + lab):
                        filename = directory + "/" + cell_type + "/" + species + "/" + region + '/' + lab + '/' + neuron
                        if neuron[-4:] != ".swc": 
                            continue
                        name = neuron[:-4]
                        outdir = 'figs/%s/%s/%s/%s' % (cell_type, species, region, lab)
                        outdir = outdir.replace(' ', '_')
                        os.system('mkdir -p %s' % outdir)

                        print species, lab, neuron
                        pareto_plot(filename, name, cell_type, species, region, lab, outdir)


if __name__ == '__main__':
    points = [(0, 0), (1, 1), (1, 1.1), (0, 0.1), (2, 2), (-1, -1), (-1, -1.1), (-1, 2), (-0.5, -0.5), (-0.5, 0.5), (0.5, 0.5), (1.1, 0.01)]
    root = (0, 0)
    #pareto_plot(points, root, 'test')

    #filename = argv[1]
    frog_filename = 'neuromorpho/frog/birinyi/GEN1.CNG.swc'
    goldfish_filename= 'neuromorpho/goldfish/stevens/G4-19g-1.CNG.swc'
    pig_filename = 'neuromorpho/pig/johnson/Pig288-DG-In-Neuron1.CNG.swc'
    agouti_filename = 'neuromorpho/agouti/guerra da rocha/cco6lam06cel05pa.CNG.swc'
    celegans_filename = 'neuromorpho/celegans/openworm/SABVR.CNG.swc' 
    
    #points = P2Coord.values()

    #print P2Coord.keys()

    #root_point = P2Coord[root]

    #pareto_plot(goldfish_filename, 'goldfish', 'cell_type1', 'species1', 'region1', 'lab1', 'figs/')
    #pareto_plot(pig_filename, 'pig')
    #pareto_plot(agouti_filename, 'agouti')
    #pareto_plot(celegans_filename, 'celegans')
    #pareto_plot(frog_filename, 'frog')
    #neuromorpho_plots()

    pareto_drawings(goldfish_filename, 'goldfish')
