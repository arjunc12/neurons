import numpy as np
from scipy.spatial.distance import euclidean
import pandas as pd

def point_dist(p1, p2):
    '''
    '''
    assert len(p1) == len(p2)
    sq_dist = 0
    for i in xrange(len(p1)):
        x1, x2 = p1[i], p2[i]
        sq_dist += (x1 - x2) ** 2
    return sq_dist ** 0.5
    '''
    '''
    #return euclidean(p1, p2)
    #return np.linalg.norm(p1 - p2)

def node_dist(G, u, v):
    p1 = G.node[u]['coord']
    p2 = G.node[v]['coord']
    return point_dist(p1, p2)

def root_dist(G, u):
    assert 'root' in G.graph
    return node_dist(G, u, G.graph['root'])

def sort_neighbors(G):
    for u in G.nodes_iter():
        G.node[u]['close_neighbors'] = sorted(G.neighbors(u), key = lambda v : G[u][v]['length'])
    G.graph['sorted'] = True

def pareto_dist_l2(pareto_mcosts, pareto_scosts, mcost, scost):
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

def pareto_dist_scale(pareto_mcosts, pareto_scosts, mcost, scost):
    assert len(pareto_mcosts) == len(pareto_scosts)
    best_scale = float("inf")
    best_index = None
    for i in xrange(len(pareto_mcosts)):
        pareto_mcost = pareto_mcosts[i]
        assert pareto_mcost >= 0
        if pareto_mcost == 0:
            continue
        pareto_scost = pareto_scosts[i]
        assert pareto_scost >= 0
        if pareto_scost == 0:
            continue
        mcost_scale = mcost / pareto_mcost
        scost_scale = scost / pareto_scost
        scale = max(mcost_scale, scost_scale)
        if scale < best_scale:
            best_scale = scale
            best_index = i
    return best_scale, best_index

def main():
    mcost, scost = 24.105986, 381.100000
    pareto_front = pd.read_csv('/iblsn/data/Arjun/neurons/pareto_steiner_output/pareto_fronts_synthetic/140826_02_recon-13/apical_dendrite/pareto_front.csv', skipinitialspace=True)
    pareto_mcosts = pareto_front['mcost']
    pareto_scosts = pareto_front['scost']

    print pareto_dist_scale(pareto_mcosts, pareto_scosts, mcost, scost)

if __name__ == '__main__':
    main()
