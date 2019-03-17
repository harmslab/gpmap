__doc__ = """Utility functions for managing genotype-phenotype map data
and conversions.

Glossary:
--------
mutations : doct
    keys are site numbers in the genotypes. Values are alphabet of mutations at
    that sites

encoding : dict
    keys are site numbers in genotype. Values are dictionaries mapping each
    mutation to its binary representation.



"""

# -------------------------------------------------------
# Miscellaneous Python functions for random task
# -------------------------------------------------------

import itertools as it
import numpy as np
from scipy.misc import comb
from collections import OrderedDict
import warnings

# -------------------------------------------------------
# Mutation alphabets
# -------------------------------------------------------

DNA = ["A", "C", "G", "T"]

AMINO_ACIDS = ["D", "T", "S", "E", "P", "G", "A", "C", "V", "M", "I",
               "L", "Y", "F", "H", "K", "R", "W", "Q", "N"]

# -------------------------------------------------------
# Wrappers for methods that use optional imports
# -------------------------------------------------------


def ipywidgets_missing(function):
    """Wrapper checks that ipython widgets are install before trying to
    render them.
    """
    def wrapper(*args, **kwargs):
        try:

            import ipywidgets
            return function(*args, **kwargs)

        except ImportError:
            warnings.filterwarnings("once")
            warnings.warn(
                """Looks like `ipywidgets` is not installed, so widgets can't "
                "be constructed. Install before using this method.""",
                ImportWarning)

    return wrapper

# -------------------------------------------------------
# Useful methods for genotype-phenotype spaces
# -------------------------------------------------------


def get_base(logbase):
    """Get base from logbase
    Parameters
    ----------
    logbase : callable
        logarithm function
    Returns
    -------
    base : float
        returns base of logarithm.
    """
    testval = 10
    return np.exp(np.log(testval) / logbase(testval))


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences """
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def sample_phenotypes(phenotypes, errors, n=1):
    """Generate `n` phenotypes from from normal distributions. """
    samples = np.random.randn(len(phenotypes), n)
    # Apply phenotype scale and variance
    for i in range(n):
        samples[:, i] = np.multiply(samples[:, i], errors) + phenotypes
    return samples

# -------------------------------------------------------
# Utilities for searching sequence space
# -------------------------------------------------------


def find_differences(s1, s2):
    """Return the index of differences between two sequences."""
    indices = list()
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            indices.append(i)
    return indices


def farthest_genotype(reference, genotypes):
    """Find the genotype in the system that differs at the most sites. """
    mutations = 0
    for genotype in genotypes:
        differs = hamming_distance(genotype, reference)
        if differs > mutations:
            mutations = int(differs)
            mutant = str(genotype)
    return mutant

# -------------------------------------------------------
# Space enumerations
# -------------------------------------------------------


def list_binary(length):
    """List all binary strings with given length.
    """
    return np.array(["".join(seq) for seq in it.product("01", repeat=length)])


def mutations_to_encoding(wildtype, mutations):
    """ Encoding map for genotype-to-binary

    Parameters
    ---------
    wildtype: str
        Wildtype sequence.
    mutations: dict
        Mapping of each site's mutation alphabet.
        {site-number: [alphabet]}

    Returns
    -------
    encode: OrderedDict of OrderDicts
        Encoding dictionary that maps site number to mutation-binary map


    Examples
    --------
    ``{ <site-number> : {<mutation>: <binary>} }``
    """
    encoding = OrderedDict()

    for site_number, alphabet in mutations.items():
        site_number = int(site_number)
        # Handle sites that don't mutate.
        if alphabet is None:
            encoding[site_number] = wildtype[site_number]

        # All sites that mutate, give a mapping dictionary.
        else:
            # copy alphabet to avoid removing items in main object.
            alphabet_cp = alphabet[:]
            n = len(alphabet_cp) - 1  # number of mutation neighbors
            wt_site = wildtype[site_number]  # wildtype letter

            # Build a binary representation of mutation alphabet
            indiv_encode = OrderedDict({wt_site: "0" * n})
            alphabet_ = alphabet_cp[:]
            alphabet_.remove(wt_site)

            for i in range(n):
                binary = list("0" * n)
                binary[i] = "1"
                indiv_encode[alphabet_[i]] = "".join(binary)
            encoding[site_number] = indiv_encode

    return encoding


def mutations_to_genotypes(mutations, wildtype=None):
    """Use a mutations dictionary to construct an array of genotypes composed
    of those mutations.

    Parameters
    ----------
    mutations : dict
        A mapping dict with site numbers as keys and lists of mutations as
        values.

    wildtype : str
        wildtype genotype (as string).

    Returns
    -------
    genotypes : list
        list of genotypes comprised of mutations in given dictionary.
    """
    # Convert mutations dict to list of lists
    mutations_ = []
    for i, val in enumerate(mutations.values()):
        if val is None:
            mutations_.append(wildtype[i])
        else:
            mutations_.append(val)
    sequences = it.product(*mutations_)
    genotypes = ["".join(s) for s in sequences]
    return genotypes


def genotypes_to_mutations(genotypes):
    """Create mutations dictionary from a list of mutations.
    """
    # Sequences to array
    arr = np.array([list(g) for g in genotypes])

    # Mutations dict
    mutations = {i: None for i in range(arr.shape[1])}

    # Find unique residues at all sites.
    for i, col in enumerate(arr.T):
        mutations[i] = list(np.unique(col))

    return mutations


def genotypes_to_binary(wildtype, genotypes, mutations):
    """Get binary representation of genotypes w.r.t. to wildtype.

    Parameters
    ----------
    wildtype : str
        wildtype sequence.

    genotypes : list
        List of genotypes to transform.

    mutations : dict
        mutations dictionary that maps sites to mutations.

    Returns
    -------
    binary : list
        list of binary representations.
    """
    # ---------- Sanity Checks ---------------
    # 1. Check genotypes are all same length
    length_of_genotypes = [len(g) for g in genotypes]
    length = length_of_genotypes[0]

    if len(set(length_of_genotypes)) > 1:
        raise Exception("Genotypes are not all the same length.")

    if len(wildtype) != length:
        raise Exception("Wildtype is not the same length as genotypes.")

    if len(mutations) != length:
        raise Exception("mutations dict is not the same length as genotypes.")

    # Encoding dictionary
    encoding = mutations_to_encoding(wildtype, mutations)

    binary = []
    for g in genotypes:
        b = ''
        for site, mutation in enumerate(g):
            if mutations[site] is not None:
                b += encoding[site][mutation]
        binary.append(b)
    return binary


def get_missing_genotypes(genotypes, mutations=None):
    """Get a list of genotypes not found in the given genotypes list.

    Parameters
    ----------
    genotypes : list
        List of genotypes.

    mutations : dict (optional)
        Mutation dictionary

    Return
    ------
    missing_genotypes : list
        List of genotypes not found in genotypes list.
    """
    if mutations is None:
        mutations = genotypes_to_mutations(genotypes)

    # Need a wildtype--doesn't matter what it is.
    wildtype = "".join([sites[0] for sites in mutations.values()])

    # Get all genotypes.
    all_genotypes = mutations_to_genotypes(mutations, wildtype)

    # Find genotypes not found in genotypes list.
    missing_genotypes = set(all_genotypes).difference(set(genotypes))
    return list(missing_genotypes)

def length_to_mutations(length, alphabet=["0", "1"]):
    """Build a mutations dictionary for a given alphabet
    
    Parameters
    ----------
    length : int
        length of the genotypes

    alphabet : list
        List of mutations at each site.  
    """
    return {i: alphabet for i in range(length)}