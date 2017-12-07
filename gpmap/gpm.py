#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Outside imports
# ----------------------------------------------------------

import json
import pickle
import numpy as np
import pandas as pd

# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

# import different maps into this module
import gpmap.mapping as mapping
import gpmap.utils as utils
import gpmap.sample as sample
import gpmap.errors as errors
import gpmap.binary as binary


class GenotypePhenotypeMap(mapping.BaseMap):
    """
    """
    def __init__(self, wildtype, genotypes, phenotypes,
                 stdeviations=None,
                 mutations=None,
                 n_replicates=1,
                 include_binary=True,
                 **kwargs):

        # Set mutations; if not given, assume binary space.
        if mutations is not None:
            # Make sure the keys in the mutations dict are integers, not
            # strings.
            self._mutations = dict([(int(key), val)
                                   for key, val in mutations.items()])
        else:
            mutant = utils.farthest_genotype(wildtype, genotypes)
            mutations = utils.binary_mutations_map(wildtype, mutant)
            self._mutations = mutations

        # Set wildtype.
        self._wildtype = wildtype

        # Store data in DataFrame
        data = dict(
            genotypes=genotypes,
            phenotype=phenotypes,
            n_replicates=n_replicates,
            stdeviations=stdeviations
        )
        self.data = pd.DataFrame(data)

        # Built the binary representation of the genotype-phenotype.
        # Constructs a complete sequence space and stores genotypes missing
        # in the data as an attribute, `missing_genotypes`.
        if include_binary:
            self.add_binary(self.wildtype)

        # Construct the error maps
        self._add_error()

    @property
    def length(self):
        """Get length of the genotypes. """
        return len(self.wildtype)

    @property
    def n(self):
        """Get number of genotypes, i.e. size of the genotype-phenotype map."""
        return len(self.genotypes)

    @property
    def wildtype(self):
        """Get reference genotypes for interactions. """
        return self._wildtype

    @property
    def mutant(self):
        """Get the farthest mutant in genotype-phenotype map."""
        _mutant = []
        _wt = self.wildtype
        for i in range(0, len(self.mutations)):
            site = _wt[i]
            options = self.mutations[i]
            if options is None:
                _mutant.append(_wt[i])
            else:
                for o in options:
                    if o != site:
                        _mutant.append(o)
        return "".join(_mutant)

    @property
    def mutations(self):
        """Get the furthest genotype from the wildtype genotype."""
        return self._mutations

    @property
    def genotypes(self):
        """Get the genotypes of the system."""
        return self.data.genotypes

    @property
    def missing_genotypes(self):
        """Genotypes that are missing from the complete genotype-to-phenotype
        map."""
        return self.missing_data.genotypes

    @property
    def complete_genotypes(self):
        """Array of sorted genotypes for the complete genotype space encoded by
        the mutations dictionary.

        **NOTE** Can only be set by the BinaryMap object.
        """
        try:
            return self.complete_data.genotypes
        except AttributeError:
            raise AttributeError("Looks like a BinaryMap has not been built "
                                 "yet for this map. Do this before asking for "
                                 "the complete_genotypes.")

    @property
    def phenotypes(self):
        """Get the phenotypes of the system. """
        return self.data.phenotypes

    @property
    def stdeviations(self):
        """Get stdeviations"""
        return self.data.stdeviations

    @property
    def n_replicates(self):
        """Return the number of replicate measurements made of the phenotype"""
        return self.data.n_replicates

    @property
    def index(self):
        """Return numpy array of genotypes position. """
        return self.data.index

    def _add_error(self):
        """Store error maps"""
        self.std = errors.StandardDeviationMap(self)
        self.err = errors.StandardErrorMap(self)

    def add_binary(self, wildtype):
        """Add a BinaryMap to the GenotypePhenotypeMap. The wildtype determines
        the encoding pattern. Wildtype sites are represented as 0's.
        """
        self.binary = binary.BinaryMap(self, wildtype)
