# Mapping object for tracking a binary representation of the epistasis map.
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------
import numpy as np

from seqspace.base import BaseMap
from seqspace.errors import StandardErrorMap, StandardDeviationMap
from seqspace.utils import (hamming_distance,
                            binary_mutations_map,
                            farthest_genotype,
                            encode_mutations,
                            construct_genotypes)

class BinaryMap(BaseMap):
    """Object holds binary representation of phenotype.

    Parameters
    ----------
    GPM : GenotypePhenotypeMap object
        The genotype phenotype map object to translate as Binary.

    Attributes
    ----------
    length : int
        length of the binary sequences
    genotypes : np.array
        array of binary genotype strings, ordered the same as input from GPM.
    missing_genotypes : np.array
        other genotypes possible by mutations, not given in the data. These are
        often the genotypes to predict.
    complete_genotypes : np.array
        genotypes + missing_genotypes
    phenotypes : np.array
        phenotypes given by GPM, in the same order as GPM.
    encoding : dict
        mapping dictionary that takes
    n_replicates : int
        number of replicates.
    logbase : callable
        function for log transforming an array or value.
    stdeviations : array
        standard deviations of genotype phenotype map.
    """
    def __init__(self, GPM):
        self._GPM = GPM
        self._build()
        self.std = StandardDeviationMap(self)
        self.err = StandardErrorMap(self)

    @property
    def n_replicates(self):
        """Get number of replicates"""
        return self._GPM.n_replicates

    @property
    def logbase(self):
        """Get logbase."""
        return self._GPM.logbase

    @property
    def stdeviations(self):
        """"""
        return self._GPM.stdeviations

    @property
    def tranformed(self):
        """Get boolean for tranformed"""
        return self._GPM.transformed

    @property
    def length(self):
        """ Get length of binary strings in space. """
        return self._length

    @property
    def genotypes(self):
        """ Get Binary representation of genotypes. """
        return self._genotypes

    @property
    def phenotypes(self):
        """Get phenotypes of the map."""
        return self._GPM.phenotypes

    @property
    def missing_genotypes(self):
        """ Binary genotypes missing in the dataset """
        return self._missing_genotypes

    @property
    def complete_genotypes(self):
        """ All possible genotypes in the complete genotype space"""
        return np.concatenate((self.genotypes, self.missing_genotypes))

    # ----------------------------------------------------------
    # Setter methods
    # ----------------------------------------------------------

    @genotypes.setter
    def genotypes(self, genotypes):
        """ Set Binary representation of genotypes. """
        self._length = len(genotypes[0])
        self._genotypes = genotypes

    def _build(self):
        """Builds a binary representation of a GenotypePhenotypeMap object.
        """
        self.encoding = encode_mutations(self._GPM.wildtype, self._GPM.mutations)

        # Use encoding map to construct binary presentation for any type of alphabet
        unsorted_genotypes, unsorted_binary = construct_genotypes(self.encoding)

        # length of binary strings
        self._length = len(unsorted_binary[0])

        # Sort binary representation to match genotypes
        mapping = self._GPM.map("genotypes", "indices")
        binary = np.empty(self._GPM.n, dtype="<U" + str(self.length))

        # Sort the genotypes by looking for them in the data.
        missing_genotypes = list()
        missing_binary = list()

        for i in range(len(unsorted_genotypes)):
            # Keep and sort genotype if it exists in data.
            try:
                index = mapping[unsorted_genotypes[i]]
                binary[index] = unsorted_binary[i]
            # If the genotype is not there.
            except KeyError:
                missing_genotypes.append(unsorted_genotypes[i])
                missing_binary.append(unsorted_binary[i])

        # Set the missing genotypes
        self._GPM._missing_genotypes = np.array(missing_genotypes)
        self._missing_genotypes = np.array(missing_binary)

        # Set binary attributes to sorted genotypes
        self._genotypes = binary
