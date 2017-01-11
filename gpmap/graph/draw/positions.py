import numpy as np

def flattened(G, scale=1, vertical=False):
    """Get flattened positions for a genotype-phenotype graph.

    Parameters
    ----------
    G : GenotypePhenotypeGraph object
        A genotype-phenotype objects
    scale : float (default=1)
        density of the nodes.

    Returns
    -------
    positions: dict
        positions of all nodes in network (i.e. {index: [x,y]})
    """
    # Get the binary genotypes from GPM
    # Set level of nodes and begin calc offset on the fly
    graph = G
    offsets = {}
    positions = {}
    for n in range(len(list(G.nodes()))):
        node = graph.node[n]
        # Calculate the level of each node
        level = node["binary"].count("1")
        if level in offsets:
            offsets[level] += 1
        else:
            offsets[level] = 1
        positions[n] = [level]
    # Center the offsets on 0
    for key, val in offsets.items():
        offsets[key] = list(np.arange(val) - (val-1)/2.0)
    # Offset positions
    if vertical:
        for n in graph.nodes():
            pos = offsets[positions[n][0]].pop(0)
            scaled = scale*pos
            positions[n].insert(0, scaled)
            positions[n][-1] *= -1
    else:
        for n in graph.nodes():
            pos = offsets[positions[n][0]].pop(0)
            scaled = scale*pos
            positions[n].append(scaled)
    return positions
