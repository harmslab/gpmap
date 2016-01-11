import numpy as np
import networkx as nx
from scipy.misc import comb
from matplotlib.cm import gray, spring 
import matplotlib.pyplot as plt


def flattened_positions(space):
    """
        Draw a flattened representation of the genotype-phenotype graph.
        
        Input:
        -----
        space: GenotypePhenotypeMap object
        
        Returns:
        -------
        positions: dict
            positions of all nodes in network (i.e. {index: [x,y]})
    """
    
    # Get the binary genotypes from GPM 
    nodes = space.Binary.genotypes
    
    # Get mapping of binary genotypes to their graph indices
    mapping = space.Binary.get_map("genotypes", "indices")

    # Build an offset dictionary as we go...
    offsets = dict([(j, 0) for j in range(space.length+1)])
    
    # Init main positions dict
    positions = {}

    for i in range(space.n):
        # count the number of mutations for horizontal axis
        x = nodes[i].count("1")
        
        # Number of vertical positions
        pascal = comb(space.length, x)
        
        # Calculate the y position
        y = -(pascal-1) / 2 + offsets[x]
        
        # Add new position
        positions[mapping[nodes[i]]] = [x,y]
        
        # Iterate offset for that index on horizontal axis
        offsets[x] += 1.0
        
    return positions


def draw_trajectories(G, trajectories, pos=None):
    pass
    
# ------------------------------------------------
# Methods for drawing trajectories on networks
# ------------------------------------------------

def edge_arrows(pos, edges):
    """ Maker a list of edge arrows. """
    arrows = list()
    for e in edges:
        arrows.append((pos[e[0]][0], pos[e[0]][1], 
                       .8*(pos[e[1]][0]-pos[e[0]][0]), 
                       .8*(pos[e[1]][1]- pos[e[0]][1]), 
                       edges[e]))
    return arrows

def edge_weight(traj):
    """ Count the number of times each edge is visited. """
    edge_counter = dict()
    for t in traj:
        sequences = t
        edges = [(sequences[i-1],sequences[i]) for i in range(1,len(sequences))]
        for e in edges:
            if e in edge_counter:
                edge_counter[e] += traj[t]
            else:
                edge_counter[e] = traj[t]
    return edge_counter

def draw_space(G, pos=None):
    """ Draw the trajectories on Graph. """
    fig = plt.figure(figsize=[7,7])
    if pos is None:
        pos = nx.spring_layout(G, iterations=150)
     
    # Draw network
    colors = list()
    for n in G.nodes():
        colors.append(G.node[n]["phenotype"])
                      
    nx.draw(G,pos, alpha=.8, 
            cmap=gray, 
            node_color=colors, 
            node_size=400, 
            arrows=False, 
            with_labels=True, 
            width=0.5,
            vmin = 0.94,
            vmax = 1.2,
           )
    return pos, fig

def draw_traj(G, traj, pos=None):
    """ Draw the trajectories on Graph. """
    fig = plt.figure(figsize=[7,7])
    if pos is None:
        pos = nx.spring_layout(G, iterations=150)
     
    # Draw network
    colors = list()
    for n in G.nodes():
        colors.append(G.node[n]["phenotype"])
                      
    nx.draw(G,pos, alpha=.8, 
            cmap=gray, 
            node_color=colors, 
            node_size=400, 
            arrows=False, 
            with_labels=True, 
            width=2,
            vmin = 0.94,
            vmax = 1.2,
           )
    
    # Draw arrows
    edges = edge_weight(traj)
    arrows = edge_arrows(pos, edges)
    for a in arrows:
        plt.arrow(a[0], a[1], a[2], a[3], alpha=0.6, width=0.005*np.log(a[4]), head_width=0.05, head_length=0.05, fc='b', ec='k')
    return pos, fig
