# -------------------------------------------------------
# Miscellaneous functions for constructing GPMs
# -------------------------------------------------------

import numpy as np
from collections import OrderedDict

# -------------------------------------------------------
# Mutation alphabets
# -------------------------------------------------------

DNA = ["A", "C", "G", "T"]

AMINO_ACIDS = ["D","T", "S", "E", "P", "G", "A", "C", "V", "M", "I"
                "L", "Y", "F", "H", "K", "R", "W", "Q", "N"]

# -------------------------------------------------------
# Useful methods for genotype-phenotype spaces
# -------------------------------------------------------

def hamming_distance(s1, s2):
    """ Return the Hamming distance between equal-length sequences """
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def find_differences(s1, s2):
    """ Return the index of differences between two sequences."""
    indices = list()
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            indices.append(i)
    return indices
    
def farthest_genotype(reference, genotypes):
    """ Find the genotype in the system that differs at the most sites. """ 
    mutations = 0
    for genotype in genotypes:
        differs = hamming_distance(genotype, reference)
        if differs > mutations:
            mutations = int(differs)
            mutant = str(genotype)
    return mutant
    
def binary_mutations_map(wildtype, mutant):
    """ Construct a site-to-binary-mutations dict between two sequences. 
    
        Args:
        ----
        wildtype: str
            wildtype sequence
        mutant: str
            mutant sequence
            
        Returns:
        -------
        mutations: dict
        
        ex.
             mutations = {
                 1: [wildtype[0], mutant[0],
                 2: [wildtype[1], mutant[1]
                 ...
             }
    """
    mutations = dict()
    for i in range(len(wildtype)):
        mutations[i+1] = [wildtype[i], mutant[i]]
    return mutations

# -------------------------------------------------------
# Space enumerations
# -------------------------------------------------------

def list_binary(length):
    """ List all binary strings with given length. """
    return np.array(["".join(seq) for seq in it.product("01", repeat=length)])

def encode_mutations(wildtype, site_alphabet):
    """ Encoding map for genotype-to-binary
    
        Args:
        ----
        wildtype: str
            Wildtype sequence.
        site_alphabet: dict
            Mapping of each site's mutation alphabet.
            {site-number: [alphabet]}
        
        Returns:
        -------
        encode: OrderedDict of OrderDicts
            Encoding dictionary that maps site number to mutation-binary map
            
            Ex:
            {
                site-number: {"mutation": "binary"},
                .
                .
                .
            }
    
    """
    encoding = OrderedDict()

    for site_number, alphabet in site_alphabet.items():
        n = len(alphabet)-1 # number of mutation neighbors
        wt_site = wildtype[site_number-1] # wildtype letter

        # Build a binary representation of mutation alphabet
        indiv_encode = OrderedDict({wt_site: "0"*n})
        alphabet_ = list(alphabet)
        alphabet_.remove(wt_site)

        for i in range(n):
            binary = list("0"*n)
            binary[i] = "1"
            indiv_encode[alphabet_[i]] = "".join(binary)
        encoding[site_number] = indiv_encode
        
    return encoding


def construct_genotypes(mutation_encoding):
    """ Constructs binary representation of genotype map given a specific alphabet
        for each site.
        
        Args:
        ----
        encode: OrderedDict of OrderDicts
            Encoding dictionary that maps site number to mutation-binary map
            
            Ex:
            {
                site-number: {"mutation": "binary"},
                .
                .
                .
            }
            
        Returns:
        -------
        genotypes: array
            Array of genotypes
        binary: array
            Array of binary represention of genotypes
    """
    
    binary = [""]
    genotypes = [""]
    for site in mutation_encoding:

        # Parameters that are needed for looping
        n_genotypes = len(genotypes)
        n_copies = len(mutation_encoding[site])
        copy_genotypes = list(genotypes)
        copy_binary = list(binary)

        # Make copies of previous sites' genotypes
        # for appending next sites binary combinations
        for i in range(n_copies-1):
            genotypes += copy_genotypes
            binary += copy_binary

        # Enumerate all possible configurations to append
        # next sites binary combinations to old
        skips = 0
        for key, val in mutation_encoding[site].items():
            for i in range(n_genotypes):
                genotypes[skips*n_genotypes + i] += key            
                binary[skips*n_genotypes + i] += val
            skips += 1

    genotypes = np.array(genotypes)
    binary = np.array(binary)
    
    return genotypes, binary