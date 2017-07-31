# Mapping object for tracking a binary representation of the epistasis map.
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------
import numpy as np
import pandas as pd

from gpmap.mapping import BaseMap
from gpmap.errors import StandardErrorMap, StandardDeviationMap
from gpmap.utils import (hamming_distance,
                            binary_mutations_map,
                            farthest_genotype,
                            encode_mutations,
                            construct_genotypes)

class BinaryMap(object):
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
    def __init__(self, GPM, wildtype):
        self._GPM = GPM
        self.wildtype = wildtype
        self.std = StandardDeviationMap(self)
        self.err = StandardErrorMap(self)

    @property
    def wildtype(self):
        """Reference genotype to define the binary representation with respect to."""
        return self._wildtype

    @wildtype.setter
    def wildtype(self, wildtype):
        """Set the reference genotype and rebuild the map."""
        self._wildtype = wildtype
        self._build()

    @property
    def n_replicates(self):
        """Get number of replicates"""
        return self._GPM.n_replicates

    @property
    def stdeviations(self):
        """Get standard deviations"""
        return self._GPM.stdeviations

    @property
    def length(self):
        """Get length of binary strings in space. """
        return self._length

    @property
    def genotypes(self):
        """Get Binary representation of genotypes. """
        return self._genotypes

    @property
    def phenotypes(self):
        """Get phenotypes of the map."""
        return self._GPM.phenotypes

    @property
    def missing_genotypes(self):
        """Binary genotypes missing in the dataset """
        return self._missing_genotypes

    @property
    def complete_genotypes(self):
        """All possible genotypes in the complete genotype space."""
        return self._complete_genotypes

    # ----------------------------------------------------------
    # Setter methods
    # ----------------------------------------------------------

    def _build(self):
        """Builds a binary representation of a GenotypePhenotypeMap object.
        """
        self.encoding = encode_mutations(self._wildtype, self._GPM.mutations)
        # Use encoding map to construct binary presentation for any type of alphabet
        unsorted_genotypes, unsorted_binary = construct_genotypes(self.encoding)

        # determine length of binary strings
        self._length = len(unsorted_binary[0])

        # Series of all possible genotypes and their binary representation
        bins = pd.Series(unsorted_binary, index=unsorted_genotypes)
        # Sort data in alphabetical order using the actual genotypes (not binary representation)
        bins = bins.sort_index()
        self._complete_genotypes = bins.reset_index(drop=True)
        self._GPM._complete_genotypes = pd.Series(bins.index)

        ## Separate observed genotypes from missing genotypes
        # Observed
        genotypes = self._GPM.genotypes
        self._genotypes = bins[genotypes].reset_index(drop=True)

        # Missing
        mask = bins.isin(genotypes)
        missing_genotypes = pd.Series(bins[~mask].index)
        missing_binary = pd.Series(bins[~mask].values)
        self._GPM._missing_genotypes = missing_genotypes
        self._missing_genotypes = missing_binary
