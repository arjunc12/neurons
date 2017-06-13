import networkx as nx



def initialize_khuller(mst):
    root = H.graph['root']
    H.node[root]['parent'] = None
    H.node[root]['droot'] = 0

def pareto_khuller(G, alpha, mst=None):
    if mst == None:
        mst = nx.minimum_spanning_tree(G, weight='length')
    H = mst.copy()
    

