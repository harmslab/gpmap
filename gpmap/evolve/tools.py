import numpy as np
from .models import wright_fisher
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
        if alphabet is not None:
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
        if alphabet is not None:
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

def transition_matrix(gpm, model=wright_fisher, N=1e8, mu=1e-8):
    """Create a transition matrix for a GenotypePhenotypeMap object.
    """
    #### Need to check gpm

    T = np.zeros((gpm.n,gpm.n), dtype=float)
    # maps
    mapping_g = gpm.map("indices", "genotypes")
    mapping_g_ = gpm.map("genotypes", "indices")
    mapping_p = gpm.map("indices", "phenotypes")

    # Build transition matrix
    for i in range(gpm.n):
        source_g = mapping_g[i]
        i_p = mapping_p[i]
        target_gs = get_neighbors(source_g, gpm.mutations)
        for target in target_gs:
            j = mapping_g_[target]
            j_p = mapping_p[j]
            x = mu * N * model(i_p, j_p, N=N)
            #print(i_p, j_p, x)
            T[i][j] = mu * N * model(i_p, j_p, N=N)

    # Populate diagonal
    for i in range(gpm.n):
        T[i,i] = 1 - sum(T[i])
    return T
