import numpy as np
from gpmap.utils import hamming_distance

def get_neighbors(genotype, mutations):
    """List all neighbors that are a single a mutation away from genotype.

    Parameters
    ----------
    genotype: str
        reference genotype.
    mutations : dict
        sites (keys) mapped to an alphabet list in genotype space (values).

    Returns
    -------
    neighbors : list
        List of neighbor genotypes
    """
    sites = list(genotype)
    neighbors = []
    for i, alphabet in mutations.items():
        # Copy alphabet to avoid over-writing
        alphabet = alphabet[:]
        alphabet.remove(sites[i])
        # Replace letters
        for a in alphabet:
            g = sites[:]
            g[i] = a
            neighbors.append("".join(g))
    return neighbors

def get_forward_neighbors(source, current, mutations):
    """List all neighbors that are a single a mutation away from genotype and
    move away from the source.

    Parameters
    ----------
    source : str
        source genotype which determines the direction to be moving away.
    current: str
        reference genotype.
    mutations : dict
        sites (keys) mapped to an alphabet list in genotype space (values).

    Returns
    -------
    neighbors : list
        List of neighbor genotypes
    """
    s_sites = list(source)
    sites = list(current)
    hd = hamming_distance(source, current)
    neighbors = []
    for i, alphabet in mutations.items():
        # Copy alphabet to avoid over-writing
        alphabet = alphabet[:]
        alphabet.remove(sites[i])
        # Replace letters
        for a in alphabet:
            g = sites[:]
            g[i] = a
            if hamming_distance(source, g) > hd:
                neighbors.append("".join(g))
    return neighbors


def monte_carlo(gpm, source, target, model, max_moves=1000, forward=False, **kwargs):
    """Use a Monte Carlo approach to sample a single trajectory between a source
    genotype and a target genotype in a genotype-phenotype map.

    The only edges accessible to a given genotype in this implementation are genotypes
    that differ by a single mutation. All other moves are ignored. The algorithm
    builds a list of neighbors on-the-fly for each step. There is no `self` probability
    considered when making a move, thus, this will NOT recapitulate stationary
    frequencies, uncover a fitness landscape, or find equilibrium states. For the sake of
    efficiency, it merely samples pathways from source to target. If you'd like
    a better sampling of the fitness landscape and its equilibrium states, try
    the monte_carlo_metropolis_criterion function.

    Parameters
    ----------
    gpm : GenotypePhenotypeMap object (or subclassed object)
        The genotype-phenotype map to sample.
    source : str
        The starting genotype for simulation.
    target : str
        The ending genotype for simulation.
    model : callable
        A callable evolutionary model function that calculates the probability
        of transitioning between two genotypes.
    max_moves : int (default=1000)
        The max number of moves to try, in case the simulation gets stuck.
    forward : bool (default True)
        Set to True to only consider forward moves. Slows down the get_neighbors
        call, but avoids longer paths.

    Keyword Arguments
    -----------------
    Keyword arguments get passed directly to the model.

    Returns
    -------
    visited : tuple
        A tuple of all genotypes visited along a trajectory.
    """
    # acquire a mapping of genotype to phenotype
    mapping_p = gpm.map("genotypes", "phenotypes")
    # Set up get_neighbors method
    args = []
    if forward is True:
        args.append(source)
        neighbors_method = get_forward_neighbors
    else:
        neighbors_method = get_neighbors
    # Begin Monte Carlo loop.
    visited = (source,)
    moves = 0
    while visited[-1] != target and moves <= max_moves:
        # Observe new genotype
        current = visited[-1]
        fitness0 = mapping_p[current]
        # Find neighbors and calculate the probability of transitioning (normalized)
        nb_args = args[:] + [current, gpm.mutations]
        neighbors = np.array(neighbors_method(*nb_args))
        fitnesses = np.array([model(fitness0, mapping_p[n], **kwargs) for n in neighbors])
        norm = fitnesses.sum()
        # Check for possible moves.
        if norm == 0:
            raise Exception ("Simulation got stuck on the '" + current + "' genotype. "
            "All neighbors are deleterious. \n"
            "Path : " + str(visited))
        # Calculate a cumulative distribution to Monte Carlo sample neighbors.
        cumulative_dist = np.array([sum(fitnesses[:i+1])/norm for i in range(len(fitnesses))])
        # Monte Carlo number to sample
        mc_number = np.random.rand()
        # Make move
        new = neighbors[cumulative_dist>=mc_number][0]
        visited += (new,)
        moves += 1
    # Check for convergence and return visited.
    if moves > max_moves:
        raise Exception("Monte Carlo exceeded max number of moves.")
    return visited

def monte_carlo_metropolis_criterion(gpm, source, target, model, max_fails=1000, **kwargs):
    """Use a Monte Carlo, Metropolis Criterion method to sample a single path
    through a genotype-phenotype map.

    The only edges accessible to a given genotype in this implementation are genotypes
    that differ by a single mutation. All other moves are ignored. The algorithm
    builds a list of neighbors on-the-fly for each step. This method chooses a sample
    at random from its neighbors and uses a Metropolis criterion to accept or
    reject the move. The output will include all moves in the simulation, including
    all 'self' moves. This is useful for sampling the fitness landscape's stationary
    frequencies.

    Parameters
    ----------
    gpm : GenotypePhenotypeMap object (or subclassed object)
        The genotype-phenotype map to sample.
    source : str
        The starting genotype for simulation.
    target : str
        The ending genotype for simulation.
    model : callable
        A callable evolutionary model function that calculates the probability
        of transitioning between two genotypes.
    max_fails : int (default=1000)
        The max number of failed moves, in case the simulation gets stuck.

    Keyword Arguments
    -----------------
    Keyword arguments get passed directly to the model.

    Returns
    -------
    visited : tuple
        A tuple of all genotypes visited along a trajectory.
    """
    # acquire a mapping of genotype to phenotype
    mapping_p = gpm.map("genotypes", "phenotypes")
    # Begin Monte Carlo loop.
    visited = (source,)
    fails = 0
    while visited[-1] != target and fails <= max_fails:
        # Observe new genotype
        current = visited[-1]
        fitness0 = mapping_p[current]
        # Find neighbors and calculate the probability of transitioning (normalized)
        neighbors = np.array(get_neighbors(current, gpm.mutations))
        # sample neighbors
        mc_choice = np.random.choice(neighbors)
        mc_fitness = model(fitness0, mapping_p[mc_choice], **kwargs)
        # Metropolis criterion
        mc_number = np.random.rand()
        if mc_number < mc_fitness:
            visited += (mc_choice,)
        else:
            visited += (current,)
            fails += 1
    # Check for convergence and return visited.
    if fails > max_fails:
        raise Exception("Monte Carlo exceeded max number of moves.")
    return visited
