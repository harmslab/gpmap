import networkx as nx
import numpy as np
from collections import OrderedDict

def path_probabilities(G, paths):
    """Calculate the probabilities of a list of paths in a network G. G must have
    a transition matrix attribute.
    """
    # Build a list of probabilities
    probabilities = list()
    # Iterate through all paths in paths-list
    transition_matrix = G.transition_matrix
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
    # Return normalized probabilities. If sum(probabilities) is zero, return
    # a vector of zeros.
    if sum(probabilities) == 0 :
        return list(probabilities)
    else:
        return list(np.array(probabilities)/sum(probabilities))


def paths_and_probabilities(G, source, target, transition_model=None, *args, **kwargs):
    """Find the most likely shortest path between a source and target.

    Parameters
    ----------
    G : GenotypePhenotypeGraph (Subclass of networkx.DiGraph)
        Networkx object constructed for genotype-phenotype map.
    source : int
        index of source node.
    target : int or list of ints
        index of target node(s).
    transition_model : callable
        function to compute transition probabilities

    Notes
    -----
    ``args`` and ``kwargs`` get passed to the transition_model function.

    Returns
    -------
    paths : list
        a list of lists representing all paths from source to target.
    probabilities : list
        a list of the probabilities for each path between source and target.
    """
    # Confirm that an evolutionary model is set in the Graph.
    if transition_model is not None:
        G.add_evolutionary_model(transition_model, *args, **kwargs)
    else:
        if G.transition_model is None:
            raise Exception("""Transition model must be set or given.""")
    # Get all paths between source and target (or multiple targets)
    try:
        paths = []
        for t in target:
            paths += list(nx.all_shortest_paths(G, source, t))
    except TypeError:
        paths = list(nx.all_shortest_paths(G, source, target))
    # Build a list of probabilities
    probabilities = list()
    # Iterate through all paths in paths-list
    transition_matrix = G.transition_matrix
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
    # Return normalized probabilities. If sum(probabilities) is zero, return
    # a vector of zeros.
    if sum(probabilities) == 0 :
        return paths, list(probabilities)
    else:
        return paths, list(np.array(probabilities)/sum(probabilities))

def flux(G, source, target, transition_model=None, *args, **kwargs):
    """Calculate the probability at each edge, i.e. the flux of probability
    through each edge.

    Parameters
    ----------
    G : GenotypePhenotypeGraph (Subclass of networkx.DiGraph)
        Networkx object constructed for genotype-phenotype map.
    source : int
        index of source node.
    target : int
        index of target node.
    transition_model : callable
        function to compute transition probabilities
    """
    paths, probs = paths_and_probabilities(G, source, target,
        transition_model=transition_model,
        *args,
        **kwargs
    )
    # convert paths to tuples
    paths = [tuple(p) for p in paths]
    # map path to probability
    traj = dict(zip(paths, probs))
    # flux mapping dictionary (edge to probability)
    flux = OrderedDict([(edge,0) for edge in G.edges()])
    for path, prob in traj.items():
        # walk through trajectory and add probabilities to each edge.
        for last_i, node in enumerate(path[1:]):
            i = path[last_i]
            j = node
            flux[(i,j)] += prob
    # Return edge probabilities as array.
    return np.array(list(flux.values()))

def ml_path(G, source, target, transition_model=None, *args, **kwargs):
    """Find the most likely shortest path between a source and target.s

    Parameters
    ----------
    G : GenotypePhenotypeGraph (Subclass of networkx.DiGraph)
        Networkx object constructed for genotype-phenotype map.
    source : int
        index of source node.
    target : int
        index of target node.
    transition_model : callable
        function to compute transition probabilities

    Notes
    -----
    ``args`` and ``kwargs`` get passed to the transition_model function.

    Returns
    -------
    path : list
        a list of each node index in the ml path.
    phenotypes : list
        a list of the phenotypes for each node in ml path.
    """
    paths, probabilities = paths_and_probabilities(G, source, target,
        transition_model=transition_model, *args, **kwargs)
    index = probabilities.index(max(probabilities))
    path = paths[index]
    phenotypes =  [G.node[p]["phenotype"] for p in path]
    return path, phenotypes


def greedy_path(G, source, target):
    """Find shortest path that yields all the best moves from source to target.

    Parameters
    ----------
    G : GenotypePhenotypeGraph (Subclass of networkx.DiGraph)
        Networkx object constructed for genotype-phenotype map.
    source : int
        index of source node.
    target : int or list of ints
        index of target node(s).

    Returns
    -------
    path : list
        a list of each node index in the greedy path.
    phenotypes : list
        a list of the phenotypes for each node in greedy path.
    """
    # Determine direction of the space w.r.t. binary encoding.
    n1 = G.node[source]["binary"].count("1")
    try:
        n2 = G.node[target]["binary"].count("1")
    except TypeError:
        n2 = G.node[target[0]]["binary"].count("1")
    if n2 > n1:
        char = "1"
    else:
        char = "0"

    # Move through space from source to target
    attempts = 0
    start = source
    path = [source]
    phenotypes = [G.node[source]["phenotype"]]
    # Move along trajectory until reaching target
    while start not in target or attempts > len(G.nodes()):
        mnode = G.node[start]["binary"].count(char)
        neighbors = G.neighbors(start)
        # iterate through neighbors and find only step forward
        for n in neighbors[:]:
            node = G.node[n]
            # number of mutations from neighbor
            mneighbor = node["binary"].count(char)
            # if backwards step, remove from list
            if mneighbor < mnode:
                neighbors.remove(n)
        # get phenotypes of leftover neighbors
        vals = [G.node[n]["phenotype"] for n in neighbors]
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
