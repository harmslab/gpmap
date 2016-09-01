import networkx as nx

def probability_of_path(paths, transition_matrix):
    """Calculate the probability of a given trajectory by calculating the
    string of conditional probabilities from a transition matrix.

    Parameters
    ----------
    paths : list of lists
        List of the indicies matching the state indices in transition matrix
    transition matrix : 2d array
        Transition matrix. elements represent the probability of transitions from
        state i --> j and vice versa.

    Returns
    -------
    probabilities : array (length = len(paths))
        Probabilities for array of paths.
    """
    if type(paths) == tuple:
        paths = [paths]
    # Build a list of probabilities
    probabilities = list()
    # Iterate through all paths in paths-list
    for p in paths:
        path_length = len(p)
        # Begin by giving this path a probability of 1.
        pi = 1
        # Iterate through edges and multiply by the
        # transition probability of that edge.
        for i in range(path_length-1):
            pi *= transition_matrix[p[i],p[i+1]]
        # Append pi to probabilities
        probabilities.append(pi)
    probabilities = np.array(probabilities)
    # Return normalized probabilities. If sum(probabilities) is zero, return
    # a vector of zeros.
    if sum(probabilities) == 0 :
        return probabilities
    else:
        return probabilities/sum(probabilities)

def greedy_path(G, source, target):
    """Find shortest path that yields all the best moves from source to target.

    Parameters
    ----------
    G : Graph
        Networkx object constructed for genotype-phenotype map.
    source : int
        index of source node.
    target : int
        index of target node.
    """
    # Determine direction of the space w.r.t. binary encoding.
    n1 = Graph.node[source]["binary"].count("1")
    n2 = Graph.node[target]["binary"].count("1")
    if n1 > n2:
        char = "0"
    else:
        char = "1"

    # Move through space from source to target
    attempts = 0
    start = source
    path = [source]
    phenotypes = [Graph.node[source]["phenotype"]]
    # Move along trajectory until reaching target
    while start != target or attempts > target:
        mnode = Graph.node[start]["binary"].count(char)
        neighbors = Graph.neighbors(start)
        # iterate through neighbors and find only step forward
        for n in neighbors[:]:
            node = Graph.node[n]
            # number of mutations from neighbor
            mneighbor = node["binary"].count(char)
            # if backwards step, remove from list
            if mneighbor < mnode:
                neighbors.remove(n)
        # get phenotypes of leftover neighbors
        vals = [Graph.node[n]["phenotype"] for n in neighbors]
        # Find maximum phenotype and add to path
        phenotypes.append(max(vals))
        # Set that node as the new starting point for loop
        start = neighbors[vals.index(phenotypes[-1])]
        # Add node to path
        path.append(start)
        attempts += 1

    # Check that we reached target
    if attempts == target:
        raise Exception("""Source to target is not possible by forward trajectory""")
    return path, phenotypes
