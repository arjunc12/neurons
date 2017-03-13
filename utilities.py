#!/usr/bin/env python

from math import sqrt
import networkx as nx
import sys,re
import random

random.seed(10301949)

#==============================================================================
#                                  ERROR CHECKING
#==============================================================================
def run_checks(G):
    """ Error checking on G. """

    # Check if there are any isolated nodes.
    if nx.number_connected_components(G) > 1: assert False
        
    # Check if truly a tree.
    if G.order() != G.size() + 1: assert False


#==============================================================================
#                                   READ PLANT 
#==============================================================================
def read_neural(filename):

    # Determine arbor type and remove from filename.
    if "_axon_stp" in filename: 
        arbor_type = "2"
        filename = filename.replace("_axon_stp","")
    elif "_basaldendrite_stp" in filename: 
        arbor_type = "3"
        filename = filename.replace("_basaldendrite_stp","")
    elif "_apicaldendrite_stp" in filename: 
        arbor_type = "4"
        filename = filename.replace("_apicaldendrite_stp","")
    else:
        assert False


    # Determine axon, basaldendrite, or apicaldendrite
    G = nx.Graph()
    P2Coord = {} # u-> (x,y)
    P2Label = {} # u-> label
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
                P2Coord[root] = (float(cols[2]),float(cols[3]))

            else:
                u,v = int(cols[6]),int(cols[0])
                assert u in G and u in P2Coord
                assert not G.has_edge(u,v)
        
                P2Coord[v] = (float(cols[2]),float(cols[3]))                    
                G.add_edge(u,v)

    for u in G:
        if   u == root       : P2Label[u] = "root"
        elif G.degree(u) == 1: P2Label[u] = "tip"
        elif G.degree(u)  > 2: P2Label[u] = "branch"
        else:                  P2Label[u] = "continue"

    # Check if there are any obvious errors in G.
    run_checks(G)

    return G,P2Coord,P2Label


#==============================================================================
#                              READ STEINER TREES
#==============================================================================
def read_steiner(filename,method):
    """ Driver. """
    if method == "SmithStar":
        return read_steiner_smithstar(filename)
    elif method == "Branch":
        return read_steiner_branch(filename)
    assert False


def read_steiner_smithstar(filename):
    """ Assumes SmithStar algorithm output format. """
    # TODO: deal with seg fault.

    G = nx.Graph()
    P2Label = {} # points to label
    P2Coord = {} # points to (x,y,z)
    first_node = True
    issue = True # seg fault or something else where we don't see steiner points.
    with open(filename) as f:
        for line in f:

            # Input points.
            # > terminal[0] = Point(90.121, 64.664, 403.97)
            if line.startswith("> terminal"):                
                match = re.search('(.+)terminal\[(\d+)\](.+)',line)
                if match:
                    u = int(match.group(2))
                    G.add_node(u)

                    point = line[line.index("Point")+6:-2]
                    x,y = map(float,point.strip().split(","))                    
                    P2Coord[u] = (x,y)
                    
                    if first_node:
                        P2Label[u] = "root"
                        first_node = False
                    else:
                        P2Label[u] = "tip" # or branch; same thing.

            # Graph edges.
            # > steiner[14] adjacencies:  22 24 23
            elif "adjacencies" in line:
                match = re.search('(.+)steiner\[(\d+)\](.+)',line)
                if match:
                    u = int(match.group(2))
                    neighbors = map(int,line.strip().split(":")[1].strip().split())
                    for v in neighbors:
                        
                        assert u in G and v in G

                        G.add_edge(u,v)

            # Steiner points.
            # > steiner[14]: 83.40931 119.71333 396.92941
            elif "steiner" in line:
                issue = False
                match = re.search('(.+)steiner\[(\d+)\](.+)',line)
                if match:
                    u = int(match.group(2))
                    G.add_node(u)

                    x,y = map(float,line.strip().split(":")[1].strip().split())
                    P2Coord[u] = (x,y)                    
                    P2Label[u] = "continue"    

    # Check if there are any obvious errors in G.
    if issue: return 0,0,0

    run_checks(G)

    return G,P2Coord,P2Label



#==============================================================================
#                             COMPUTE LENGTHS
#==============================================================================
def euclidean_dist(p1,p2,P2Coord):
    """ Computes the Euclidean distance between the two points. """

    if len(P2Coord[p1]) == 2: # 2D
        x1,y1 = P2Coord[p1]
        x2,y2 = P2Coord[p2]
        
        return sqrt( (x1-x2)**2 + (y1-y2)**2 )

    elif len(P2Coord[p1]) == 3: # 3D
        x1,y1,z1 = P2Coord[p1]
        x2,y2,z2 = P2Coord[p2]

        return sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    else:
        assert False


def compute_length(G,P2Coord):
    """ Computes the total length of the plant. """
    length = sum([euclidean_dist(u,v,P2Coord) for u,v in G.edges_iter()])
        
    return length

def compute_droot(G,P2Coord,P2Label):
    """ Computes the total distance from each leaf to the root. """

    # Find the root.
    for u in G: 
        if P2Label[u] == "root":
            root = u
            break

    # Compute distance from each leaf to the root.
    droot = 0
    for u in G:
        if P2Label[u] == "tip" or P2Label[u] == "branch":
            path = nx.shortest_path(G,root,u)
            for i in xrange(len(path)-1):                
                droot += euclidean_dist(path[i],path[i+1],P2Coord)

    return droot


#==============================================================================
#                             SATELLITE ALGORITHM
#==============================================================================
def compute_sat_length(P2Coord,P2Label):
    """ The total length and total droot for the satellite solution is the same. """

    # Find the root.
    for u in P2Coord:
        if P2Label[u] == "root":
            root = u
            break

    # Compute the distance from each leaf to the root.
    length = 0
    for u in P2Coord:
        if P2Label[u] == "tip" or P2Label[u] == "branch":
            length += euclidean_dist(root,u,P2Coord)

    return length


#==============================================================================
#                               MST ALGORITHM
#==============================================================================
def compute_mst(P2Coord,P2Label):
    G = nx.Graph()
    for u in P2Coord:
        if P2Label[u] not in ["root","tip","branch"]: continue
        for v in P2Coord:
            if P2Label[v] not in ["root","tip","branch"]: continue
            if u >= v: continue
            G.add_edge(u,v,weight=euclidean_dist(u,v,P2Coord))

    return nx.minimum_spanning_tree(G)


def compute_mst_length(P2Coord,P2Label):
    
    Gmst = compute_mst(P2Coord,P2Label)
    return compute_length(Gmst,P2Coord)


def compute_mst_droot(P2Coord,P2Label):
    
    Gmst = compute_mst(P2Coord,P2Label)
    return compute_droot(Gmst,P2Coord,P2Label)



#==============================================================================
#                               RAND ALGORITHM
#==============================================================================
def compute_rand_tree(P2Coord,P2Label):
    """ Computes the lengths and droots for a random tree. """

    lengths,droots = [],[]
    for i in xrange(100):

        # Create empty graph.
        G = nx.Graph()
        for u in P2Coord:
            if P2Label[u] not in ["root","tip","branch"]: continue
            G.add_node(u)

        Isolated = G.nodes()
        Added = []
        random.shuffle(Isolated)

        # Add first edge.
        u = Isolated.pop()
        v = Isolated.pop()
        assert u != v

        G.add_edge(u,v)
        Added.append(u)
        Added.append(v)

        # Add remaining nodes.
        while len(Isolated) > 0:
            u = Isolated.pop()
            assert u not in Added
            
            v = random.choice(Added)
            assert not G.has_edge(u,v)
            
            G.add_edge(u,v)
            Added.append(v)

        run_checks(G)

        lengths.append(compute_length(G,P2Coord))
        droots.append(compute_droot(G,P2Coord,P2Label))

    return lengths,droots





    



