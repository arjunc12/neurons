import numpy as np
from scipy.spatial.distance import euclidean
from itertools import product

def closest_midpoint(P0, P1, Q0):
    P0 = np.array(P0)
    P1 = np.array(P1)
    Q0 = np.array(Q0)
    
    Pslope = P1 - P0
    a = float(np.dot(Pslope, Q0))
    b = float(np.dot(P0, Pslope))
    c = float(np.dot(Pslope, Pslope))
    t = (a - b) / c
    t = np.clip(t, 0, 1)
    return P0 + t * Pslope

def line_seg_dist(P0, P1, Q0, Q1):
    dim = len(P0)
    assert len(P1) == dim
    assert len(Q0) == dim
    assert len(Q1) == dim
    
    P0 = np.array(P0)
    P1 = np.array(P1)
    Q0 = np.array(Q0)
    Q1 = np.array(Q1)
    
    Pslope = P1 - P0
    Qslope = Q1 - Q0
    for i in xrange(dim):
        Pslope[i] = P1[i] - P0[i]
        Qslope[i] = Q1[i] - Q0[i]
            
    w0 = P0 - Q0
    
    a = float(np.dot(Pslope, Pslope))
    b = float(np.dot(Pslope, Qslope))
    c = float(np.dot(Qslope, Qslope))
    d = float(np.dot(Pslope, w0))
    e = float(np.dot(Qslope, w0))
    
    deltaPinf = None
    deltaQinf = None
    
    denom = a * c - b * b
    if denom == 0:
        deltaPinf = 0
        deltaQinf = d / float(b)
    else:
        deltaPinf = (b * e - c * d) / denom
        deltaQinf = (a * e - b * d) / denom
        
    assert deltaPinf != None
    assert deltaQinf != None
    
    candidates = []
    if (0 <= deltaPinf <= 1) and (0 <= deltaQinf <= 1):
        candidates.append((deltaPinf, deltaQinf))
    if deltaPinf < 0:
        candidates.append((0, e / c))
    elif deltaPinf > 1:
        candidates.append((1, (e + b) / c))
        
    if deltaQinf < 0:
        candidates.append((-d / a, 0))
    elif deltaQinf > 1:
        candidates.append(((b - d) / a, 1))
        
    best_dist = float("inf")
    for candidate in candidates:
        candidate = np.clip(candidate, 0, 1)
        deltaP0, deltaQ0 = candidate
        Pc = P0 + Pslope * deltaP0
        Qc = Q0 + Qslope * deltaQ0
        dist = point_dist(Pc, Qc)
        if dist < best_dist:
            best_dist = dist
            
    return best_dist

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

def k_nearest_neighbors(G, u, k=None, candidate_nodes=None):
    if candidate_nodes == None:
        candidate_nodes = list(G.nodes())
        candidate_nodes.remove(u)
    nearest_neighbors = sorted(candidate_nodes, key = lambda v : node_dist(G, u, v))
    if k != None:
        assert type(k) == int
        nearest_neighbors = nearest_neighbors[:k]
    return nearest_neighbors

def sort_neighbors(G, k=None):
    for u in G.nodes():
        #G.node[u]['close_neighbors'] = sorted(G.neighbors(u), key = lambda v : G[u][v]['length'])
        #G.node[u]['close_neighbors'] = sorted(G.nodes(), key=lambda v : node_dist(G, u, v))
        #G.node[u]['close_neighbors'].remove(u)
        G.node[u]['close_neighbors'] = k_nearest_neighbors(G, u, k=k, candidate_nodes=None)
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
    P0 = (0, 0, 0)
    P1 = (0, 0, 0.1)
    Q0 = (0, 1, 0)
    Q1 = (1, 1, 0)
    print line_seg_dist(P0, P1, Q0, Q1)
    
if __name__ == '__main__':
    main()
