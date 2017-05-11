from pareto_mst import centroid, viz_tree
import networkx as nx

G = nx.Graph()
for u in [(1, 1), (1, 2), (2, 1), (2, 2), (0, 0)]:
    G.add_node(u)
    G.node[u]['coord'] = u

G.graph['root'] = (0, 0)
print centroid(G)
