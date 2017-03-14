#!/usr/bin/env python

from __future__ import division
import matplotlib as mpl
mpl.use('agg')
from matplotlib import pylab as PP
import networkx as nx
import utilities as UTIL
import os,sys

""" Plots Total MST length vs D(root) for Masland RGC cells. """

#==============================================================================
#                                VISUALIZATION
#==============================================================================
def viz_tree(G,P2Coord,P2Label):
    """ Displays plant/tree visualization. """
    
    node_size,node_color = [],[]
    pos = {}    
    for u in G:
        
        pos[u] = (P2Coord[u][0],P2Coord[u][1])
        
        if P2Label[u] == "root":
            node_color.append('black')
            node_size.append(350)            
        elif P2Label[u] == "branch" or P2Label[u] == "tip":
            node_color.append('green')
            node_size.append(100)
        elif P2Label[u] == "continue":
            node_color.append('brown')
            #node_size.append(250)   
            node_size.append(1)   
        else:
            assert False

    nx.draw(G,pos=pos,arrows=False,with_labels=False,node_size=node_size,node_color=node_color,edge_color="brown",width=4,font_size=12,font_color='red',font_weight='bold')
    #nx.draw(G,pos=pos,arrows=False,with_labels=False,node_size=15,font_size=12,font_color='red',font_weight='bold')
    PP.draw()
    PP.show()
    #PP.savefig("qd15.pdf")
    PP.close()


#==============================================================================
#                                 DO IT!
#==============================================================================
def doit(filename,annot,dim="3D"):

    for arbor_type in ["2","3","4"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
    #for arbor_type in ["3"]: # 2 = axon, 3 = basal dendrite, 4 = apical dendrite.
        
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
                    
                    if dim == "3D":                
                        P2Coord[root] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        P2Coord[root] = (float(cols[2]),float(cols[3]))
                    else:
                        assert False


                else:
                    u,v = int(cols[6]),int(cols[0])
                    assert u in G and u in P2Coord
                    assert not G.has_edge(u,v)
            
                    if dim == "3D":
                        P2Coord[v] = (float(cols[2]),float(cols[3]),float(cols[4]))
                    elif dim == "2D":
                        P2Coord[v] = (float(cols[2]),float(cols[3]))
                    else: 
                        assert False

                    G.add_edge(u,v)

                    #if UTIL.euclidean_dist(u,v,P2Coord) > 10:
                    #    print u,v
        
        
        # Populate labels, roots, branches, and tip points.
        num_pts = 0
        P2Label = {} # u-> label        
        for u in G:
            if u == root: 
                P2Label[u] = "root"
                num_pts += 1
            elif G.degree(u) == 1: 
                P2Label[u] = "tip"
                num_pts += 1
            elif G.degree(u) > 2: 
                P2Label[u] = "branch"
                num_pts += 1
            else:                  
                P2Label[u] = "continue"

        if num_pts <= 4 or num_pts >= 1000:
            return

        UTIL.run_checks(G)

        #return G,P2Coord,P2Label
        
        # Compute neural length and d(root).
        neural_length = UTIL.compute_length(G,P2Coord)
        neural_droot  = UTIL.compute_droot(G,P2Coord,P2Label)
        
        # Compute MST length and d(root).
        mst_length = UTIL.compute_mst_length(P2Coord,P2Label)
        mst_droot  = UTIL.compute_mst_droot(P2Coord,P2Label)

        sat_length = UTIL.compute_sat_length(P2Coord,P2Label)
        sat_droot  = sat_length

        # Check if the expected Pareto relation holds.
        check = (mst_length < neural_length < sat_length) and (sat_droot < neural_droot < mst_droot)

        # Tesht versus random tree.
        rand_length,rand_droot = UTIL.compute_rand_tree(P2Coord,P2Label)
        inbox,beatneural = 0,0
        for L,D in zip(rand_length,rand_droot):
            
            beatneural += (L < neural_length) or (D < neural_droot) # beats neural in > 0 dimensions.
            inbox      += (mst_length < L < sat_length) and (sat_droot < D < mst_droot) # in pareto box.           
            
        beatneural = beatneural / len(rand_length)
        inbox = inbox / len(rand_length)
        
        # Compute ratio from opt for the plant.
        neural_ab = (neural_length / mst_length) / (neural_droot / sat_droot)
        mst_ab    = (mst_length / mst_length) / (mst_droot / sat_droot)
        sat_ab    = (sat_length / mst_length) / (sat_droot / sat_droot)

        
        # Output results.
        print "%s\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%i\t%.3f\t%.3f" %(check,neural_ab,mst_ab,sat_ab,neural_length,mst_length,sat_length,neural_droot,mst_droot,sat_droot,annot,num_pts,inbox,beatneural)#[filename.rindex("/")+1:])
        return neural_length, neural_droot, mst_length, sat_droot

def main(): 
    # OLD:
    # Iterate through masland neurons.
    # directory = "../data/neurons/neuromorpho/masland-mouse-rgc"    
    # for filename in os.listdir(directory):
    #     if filename.startswith("."): continue
    #     doit(directory+"/"+filename,"masland-"+filename,dim="3D")

    # return

    neural_dists = []
    neural_droots = []
    mst_lengths = []
    sat_droots = []
    ratios = []

    if sys.argv[1] == "rgc":

        directory = "../data/neurons/rgc"
        for lab in os.listdir(directory):
            if lab.startswith("."): continue
            for neuron in os.listdir(directory + "/" + lab):
                if neuron.startswith("."): continue

                doit(directory+"/"+lab+"/"+neuron,lab+"-"+neuron,dim="3D")
                #return

    elif sys.argv[1] == "all2d":

        # Iterate through every species, lab, and neuron.
        #directory = "../data/neurons/neuromorpho"
        directory = 'neuromorpho'
        i = 0
        for species in os.listdir(directory):
            for lab in os.listdir(directory + "/" + species):
                for neuron in os.listdir(directory + "/" + species + "/" + lab):
                    filename = directory + "/" + species + "/" + lab + "/" + neuron
                    if neuron[-4:] != ".swc": continue

                    result = doit(filename,species+"-"+lab+"-"+neuron,dim="2D")
                    if result != None:
                        i += 1
                        neural_dist, neural_droot, mst_length, sat_droot = result
                        neural_dists.append(neural_dist)
                        neural_droots.append(neural_droot)
                        mst_lengths.append(mst_length)
                        sat_droots.append(sat_droot)
                        ratios.append(float(neural_dist) / float(neural_droot))
    else:
        assert False


    PP.scatter(neural_dists, neural_droots, c='r')
    PP.scatter(mst_lengths, sat_droots, c='b')
    PP.scatter(range(len(ratios)), sorted(ratios), c='g')
    PP.savefig('neural_trees.pdf', format='pdf')
    PP.close()


if __name__ == "__main__":
    main()
