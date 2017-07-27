import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import networkx as nx

def _expand(G, explored_nodes, explored_edges):
    """
    Expand existing solution by a process akin to BFS.

    Arguments:
    ----------
    G: networkx.Graph() instance
        full graph

    explored_nodes: set of ints
        nodes visited

    explored_edges: set of 2-tuples
        edges visited

    Returns:
    --------
    solutions: list, where each entry in turns contains two sets corresponding to explored_nodes and explored_edges
        all possible expansions of explored_nodes and explored_edges

    """
    frontier_nodes = list()
    frontier_edges = list()
    for v in explored_nodes:
        for u in nx.neighbors(G,v):
            if not (u in explored_nodes):
                frontier_nodes.append(u)
                frontier_edges.append([(u,v), (v,u)])

    return zip([explored_nodes | frozenset([v]) for v in frontier_nodes], [explored_edges | frozenset(e) for e in frontier_edges])

def find_all_spanning_trees(G, root=0):
    """
    Find all spanning trees of a Graph.

    Arguments:
    ----------
    G: networkx.Graph() instance
        full graph

    Returns:
    ST: list of networkx.Graph() instances
        list of all spanning trees

    """

    # initialise solution
    explored_nodes = frozenset([root])
    explored_edges = frozenset([])
    solutions = [(explored_nodes, explored_edges)]
    # we need to expand solutions number_of_nodes-1 times
    for ii in range(G.number_of_nodes()-1):
        # get all new solutions
        solutions = [_expand(G, nodes, edges) for (nodes, edges) in solutions]
        # flatten nested structure and get unique expansions
        solutions = set([item for sublist in solutions for item in sublist])

    return [nx.from_edgelist(edges) for (nodes, edges) in solutions]


if __name__ == "__main__":

    N = 3
    G = nx.grid_2d_graph(N,N)
    labels = dict( ((i,j), i + (N-1-j) * N ) for i, j in G.nodes() )
    nx.relabel_nodes(G,labels,False)
    inds=labels.keys()
    vals=labels.values()
    inds=[(N-j-1,N-i-1) for i,j in inds]
    pos2=dict(zip(vals,inds))

    #fig, ax = plt.subplots(1,1)
    #nx.draw_networkx(G, pos=pos2, with_labels=True, node_size = 200, node_color='orange',font_size=10,ax=ax)
    #plt.axis('off')
    #plt.title('grid')

    ST = find_all_spanning_trees(G)
    print len(ST)
    
    '''
    for g in ST:
        fig, ax = plt.subplots(1,1)
        nx.draw_networkx(g, pos=pos2, with_labels=True, node_size = 200, node_color='orange',font_size=10,ax=ax)
        plt.axis('off')
        plt.title('grid')
        plt.show()
    '''
