from scipy.misc import comb

def flattened_positions(graph, scale=1, vertical=False):
    """Get flattened positions for a genotype-phenotype graph.

    Parameters
    ----------
    graph : GenotypePhenotypeGraph object
        A genotype-phenotype objects

    Returns
    -------
    positions: dict
        positions of all nodes in network (i.e. {index: [x,y]})
    """
    # Get the binary genotypes from GPM
    nodes = graph.gpm.binary.genotypes
    # Get mapping of binary genotypes to their graph indices
    mapping = graph.gpm.map("binary.genotypes", "indices")
    # Build an offset dictionary as we go...
    offsets = dict([(j, 0) for j in range(graph.gpm.length+1)])
    # Init main positions dict
    positions = {}
    for i in range(graph.gpm.n):
        # count the number of mutations for horizontal axis
        x = nodes[i].count("1")
        # Number of vertical positions
        pascal = comb(graph.gpm.length, x)
        # Calculate the y position
        y = -(pascal-1) / 2 + offsets[x]
        # Add new position
        if vertical:
            positions[mapping[nodes[i]]] = [y, -x*scale]
        else:
            positions[mapping[nodes[i]]] = [x*scale,y]
        # Iterate offset for that index on horizontal axis
        offsets[x] += 1.0
    return positions
