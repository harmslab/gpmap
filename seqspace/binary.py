# Mapping object for tracking a binary representation of the epistasis map.
#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------
import numpy as np

from seqspace.base import BaseMap
from seqspace.errors import ErrorMap

class BinaryMap(BaseMap):
    """
        Map for holding a binary representation of an epistasis map.
    """
    
    def __init__(self):
        self._errors = ErrorMap()

    @property
    def length(self):
        """ Get length of binary strings in space. """
        return self._length

    @property
    def genotypes(self):
        """ Get Binary representation of genotypes. """
        return self._genotypes

    @property
    def missing_genotypes(self):
        """ Binary genotypes missing in the dataset """
        return self._missing_genotypes

    @property
    def complete_genotypes(self):
        """ All possible genotypes in the complete genotype space"""
        return np.concatenate((self.genotypes, self.missing_genotypes))

    @property
    def indices(self):
        """ Get indices of genotypes in self.genotypes that mapped to their binary representation.

            **NOTE** This will probably change -- these indices describe how the non-binary
            arrays would need to be rearrange to fit the binary array. It would make MORE sense
            to flip this around. These indices should describe where this
            binary representation is in the non-binary array (i.e. how should the binary array be
            rearranged to align with non-binary). Further, these should be aligned from the start!

        """
        return self._indices

    @property
    def phenotypes(self):
        """ Get the phenotype values in an array orderd same as binary reprentation."""
        return self._phenotypes

    @property
    def encoding(self):
        """ Return a binary representation of each site-mutation in the genotype-phenotype map"""
        return self._encoding

    @property
    def errors(self):
        """ Get the phenotype values in an array orderd same as binary reprentation."""
        return self._errors

    # ----------------------------------------------------------
    # Setter methods
    # ----------------------------------------------------------

    @genotypes.setter
    def genotypes(self, genotypes):
        """ Set Binary representation of genotypes. """
        self._length = len(genotypes[0])
        self._genotypes = genotypes

    @indices.setter
    def indices(self, indices):
        """ Set indices of genotypes in self.genotypes that mapped to their binary representation. """
        self._indices = indices

    @encoding.setter
    def encoding(self, encoding):
        """ Set the mapping for site-to-mutation-to-binary-representation."""
        self._encoding = encoding

    @phenotypes.setter
    def phenotypes(self, phenotypes):
        """ Set the phenotype values in an array orderd same as binary reprentation."""
        self._phenotypes = phenotypes
