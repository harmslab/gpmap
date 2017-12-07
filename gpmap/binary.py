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
    """
    """
    def __init__(self, gpm, wildtype):
        self.gpm = gpm
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
        return self.gpm.data.n_replicates

    @property
    def stdeviations(self):
        """Get standard deviations"""
        return self.gpm.data.stdeviations

    @property
    def length(self):
        """Get length of binary strings in space. """
        return len(self.gpm.data.binary[0])

    @property
    def genotypes(self):
        """Get Binary representation of genotypes. """
        return self.gpm.data.binary

    @property
    def phenotypes(self):
        """Get phenotypes of the map."""
        return self.gpm.phenotypes

    @property
    def missing_genotypes(self):
        """Binary genotypes missing in the dataset """
        return self.gpm.missing_data.binary

    @property
    def complete_genotypes(self):
        """All possible genotypes in the complete genotype space.

        Sorted in alphabetical order according to the
        GenotypePhenotypeMap.complete_genotypes attribute.
        """
        return self.gpm.complete_data.binary

    def _build(self):
        """Builds a binary representation of the genotypes in
        GenotypePhenotypeMap object. Also enumerates genotypes not seen in
        the genotype-phenotype map and exposes two new attributes,
        ``missing_genotypes`` and ``complete_genotypes``.

        **NOTE**: the ``complete_genotypes`` are sorted in alphabetic order.
        """
        self.encoding = encode_mutations(self._wildtype, self.gpm.mutations)
        # Use encoding map to construct binary presentation for any type of
        # alphabet
        unsorted_genotypes, unsorted_binary = construct_genotypes(
            self.encoding)

        data = {'genotypes': unsorted_genotypes,
                'binary': unsorted_binary}

        gpm = self.gpm
        gpm.complete_data = pd.DataFrame(data)
        gpm.complete_data.sort_values('genotypes', inplace=True)
        gpm.complete_data.reset_index(drop=True, inplace=True)

        # Mapping genotypes to index
        mapping = dict(zip(gpm.complete_data.genotypes,
                           gpm.complete_data.index))

        # Get index
        observed_index = [mapping[g] for g in gpm.data.genotypes]
        missing_index = set(gpm.complete_data.index).difference(observed_index)

        # Reset index of main data.
        gpm.data.index = observed_index
        # Add a column for binary representation of genotypes.
        gpm.data['binary'] = pd.Series(gpm.complete_data.binary,
                                       index=observed_index)

        # Create a dataframe for the missing data.
        gpm.missing_data = pd.DataFrame(gpm.complete_data,
                                        index=missing_index)
