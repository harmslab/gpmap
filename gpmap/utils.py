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


def binary_mutations_map(wildtype, mutant):
    """Construct a site-to-binary-mutations dict between two sequences.

    Parameters
    ----------
    wildtype: str
        wildtype sequence
    mutant: str
        mutant sequence

    Returns
    -------
    mutations: dict

    Example
    -------
    Mutations dictionary::

         mutations = {
             0: [wildtype[0], mutant[0]],
             1: [wildtype[1], mutant[1]]
             ...
         }
    """
    mutations = dict()
    for i in range(len(wildtype)):
        if wildtype[i] == mutant[i]:
            mutations[i] = None
        else:
            mutations[i] = [wildtype[i], mutant[i]]
    return mutations

# -------------------------------------------------------
# Space enumerations
# -------------------------------------------------------


def list_binary(length):
    """List all binary strings with given length.
    """
    return np.array(["".join(seq) for seq in it.product("01", repeat=length)])


def enumerate_space(wildtype, mutant, binary=True):
    """Enumerate a binary genotype list between two sequences.

    Parameters
    ----------
    wildtype: str
        Wildtype sequence as starting reference point.
    mutant: str
        Mutant sequence.

    Returns
    -------
    genotypes : list
        List of all sequence combinations between the two sequences.
    binary : list (optional)
        List of binary representations.

    Example
    -------
    if wildtype == 'AAA' and mutant == 'TTT':
        sequence space =    ['AAA','AAV','AVA','VAA','AVV','VAV','VVA','VVV']
    """

    # Check that wildtype and mutant are the same length
    if len(wildtype) != len(mutant):
        raise IndexError("ancestor_sequence and derived sequence must be the "
                         "same length.")

    # Count mutations and keep indices
    mutations = find_differences(wildtype, mutant)
    n_mut = len(mutations)
    binary_wt = "".zfill(n_mut)
    size = 2**n_mut
    rev_mutations = [mutations[i] for i in range(n_mut - 1, -1, -1)]
    mutation_map = dict(zip(mutations, range(n_mut)))

    # Enumerate mutations flipping combinations
    combinations = np.array([list(j) for i in range(1, n_mut + 1)
                             for j in it.combinations(rev_mutations, i)])
    # Initialize empty arrays
    genotypes = np.empty(size, dtype="<U" + str(len(wildtype)))
    binaries = np.empty(size, dtype="<U" + str(n_mut))
    # Population first element with wildtypes
    genotypes[0] = wildtype
    binaries[0] = binary_wt
    # Iterate through mutations combinations and build binary representations
    counter = 1
    for c in combinations:
        sequence = list(wildtype)
        b = list(binary_wt)
        for el in c:
            sequence[el] = mutant[el]   # Sequence version of mutant
            b[mutation_map[el]] = '1'            # Binary version of mutant
        genotypes[counter] = "".join(sequence)
        binaries[counter] = "".join(b)
        counter += 1

    if binary:
        return genotypes, binaries
    else:
        return genotypes


def encode_mutations(wildtype, mutations):
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


def construct_genotypes(encoding):
    """Constructs binary representation of genotype map given a specific
    alphabet for each site.

    Parameters
    ----------
    encoding: OrderedDict of OrderDicts
        Encoding dictionary that maps site number to mutation-binary map.
        *NOTE* If site does not mutate, value is set to wildtype site string
        (not dictionary).

    Returns
    -------
    genotypes: list
        list of genotypes
    binary: list
        list of binary represention of genotypes
    """
    binary = [""]
    genotypes = [""]
    for site in encoding:
        # If the site is not mutating, just add wildtype encoding
        if type(encoding[site]) is str:
            # Parameters that are needed for looping
            n_genotypes = len(genotypes)

            # Enumerate all possible configurations to append
            # wildtype site to genotypes. Binary sequences stay the same.
            for i in range(n_genotypes):
                genotypes[i] += encoding[site]

        # else, apply binary encoding scheme for that mutation.
        else:
            # Parameters that are needed for looping
            n_genotypes = len(genotypes)
            n_copies = len(encoding[site])
            copy_genotypes = list(genotypes)
            copy_binary = list(binary)

            # Make copies of previous sites' genotypes
            # for appending next sites binary combinations
            for i in range(n_copies - 1):
                genotypes += copy_genotypes
                binary += copy_binary

            # Enumerate all possible configurations to append
            # next sites binary combinations to old
            skips = 0
            for key, val in encoding[site].items():
                for i in range(n_genotypes):
                    genotypes[skips * n_genotypes + i] += key
                    binary[skips * n_genotypes + i] += val
                skips += 1

    return genotypes, binary


def mutations_to_genotypes(wildtype, mutations):
    """Use a mutations dictionary to construct an array of genotypes composed
    of those mutations.

    Parameters
    ----------
    wildtype : str
        wildtype genotype (as string).
    mutations : dict
        A mapping dict with site numbers as keys and lists of mutations as
        values.

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
    """Get binary representation of genotypes w.r.t. to wildtype."""
    # Encoding dictionary
    encoding = encode_mutations(wildtype, mutations)

    binary = []
    for g in genotypes:
        b = ''
        for site, mutation in enumerate(g):
            b += encoding[site][mutation]
        binary.append(b)
    return binary
