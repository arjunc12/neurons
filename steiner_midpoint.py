import numpy as np
from sys import argv
from itertools import combinations
from dist_functions import *

def slope_vector(p1, p2):
    assert len(p1) == len(p2)
    '''
    '''
    slope_vec = []
    for i in xrange(len(p1)):
        slope_vec.append(p2[i] - p1[i])
    return slope_vec
    '''
    '''
    #return p2 - p1

def delta_point(p1, slope, t):
    '''
    '''
    p2 = []
    assert len(p1) == len(slope)
    for i in xrange(len(p1)):
        p2.append(p1[i] + t * slope[i])
    return p2
    '''
    '''
    #return p1 + slope * t

def midpoint_cost(p2, p3, midpoint, alpha):
    a = point_dist(p3, midpoint)
    b = point_dist(midpoint, p2)

    mcost = a
    scost = a + b
    cost = pareto_cost(mcost=mcost, scost=scost, alpha=alpha)

    return cost

def best_midpoint(p1, p2, p3, alpha):
    pass

def best_midpoint_approx(p1, p2, p3, alpha):
    '''
    print "p1"
    print map(type, p1)
    print "p2"
    print map(type, p2)
    print "p3"
    print map(type, p3)
    '''

    delta = 0.01
    slope = slope_vector(p1, p2)

    best_cost = float("inf")
    best_midpoint = None
    #for dt in np.arange(delta, 1, delta):
    dt = delta
    while dt < 1: 
        #p4 = delta_point(p1, slope, dt)
        p4 = delta_point(p1, slope, float(dt))
        
        cost = midpoint_cost(p2, p3, p4, alpha)
        if cost < best_cost:
            best_cost = cost
            best_midpoint = p4

        dt += delta

    choice = 3

    cost1 = midpoint_cost(p2, p3, p1, alpha)
    if cost1 < best_cost:
        best_cost = cost1
        best_midpoint = p1
        choice = 1

    cost2 = midpoint_cost(p2, p3, p2, alpha)
    if cost2 < best_cost:
        best_cost = cost2
        best_midpoint = p2
        choice = 2

    return best_midpoint, choice

def steiner_points(p1, p2, npoints=10):
    slope = slope_vector(p1, p2)
    delta = 1.0 / (npoints + 1)
    dt = delta
    midpoints = []
    while dt < 1:
        p3 = delta_point(p1, slope, dt)
        midpoints.append(tuple(p3))
        dt += delta

    return midpoints

def colinearity_error(p1, p2, p3):
    distances = []
    points = [p1, p2, p3]
    for point1, point2 in combinations(points, 2):
        dist = point_dist(point1, point2)
        distances.append(dist)
    assert len(distances) == 3
    distances = sorted(distances)
    return distances[2] - distances[0] - distances[1]

def main():
    p1 = (0, 0, 0)
    p2 = (1, 1, 1)
    print steiner_points(p1, p2)

if __name__ == '__main__':
    main()
