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
    """Constructs a binary representation of the genotype-phenotype map. Useful
    for building networks, constructing epistasis models, and filling in
    genotype-phenotype maps.

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
        """Reference genotype."""
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
        """All possible genotypes in the complete genotype space.

        Sorted in alphabetical order according to the
        GenotypePhenotypeMap.complete_genotypes attribute.
        """
        return self._complete_genotypes

    # ----------------------------------------------------------
    # Setter methods
    # ----------------------------------------------------------

    def _build(self):
        """Builds a binary representation of the genotypes in
        GenotypePhenotypeMap object. Also enumerates genotypes not seen in
        the genotype-phenotype map and exposes two new attributes,
        ``missing_genotypes`` and ``complete_genotypes``.

        **NOTE**: the ``complete_genotypes`` are sorted in alphabetic order.
        """
        self.encoding = encode_mutations(self._wildtype, self._GPM.mutations)
        # Use encoding map to construct binary presentation for any type of
        # alphabet
        unsorted_genotypes, unsorted_binary = construct_genotypes(
            self.encoding)

        # determine length of binary strings
        self._length = len(unsorted_binary[0])

        # Series of all possible genotypes and their binary representation
        bins = pd.Series(unsorted_binary, index=unsorted_genotypes)
        # Sort data in alphabetical order using the actual genotypes
        # (not binary representation)
        bins = bins.sort_index()

        # Aliases for (some) clarity below.
        binary = self
        true = self._GPM

        # Build complete genotype map
        binary._complete_genotypes = bins.reset_index(drop=True)
        true._complete_genotypes = pd.Series(bins.index)
        mapping = {g: i for i, g in true._complete_genotypes.iteritems()}

        # Build observed genotype map
        obs_index = np.array([mapping[g] for g in true._genotypes])
        true._genotypes = pd.Series(true._complete_genotypes, index=obs_index)
        binary._genotypes = pd.Series(binary._complete_genotypes,
                                      index=obs_index)

        # Missing
        arr = np.arange(len(bins))
        missing_index = np.delete(arr, obs_index)
        true._missing_genotypes = pd.Series(true._complete_genotypes,
                                            index=missing_index)
        binary._missing_genotypes = pd.Series(binary._complete_genotypes,
                                              index=missing_index)
