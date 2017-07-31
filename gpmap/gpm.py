#
# Author: Zach Sailer
#
# ----------------------------------------------------------
# Outside imports
# ----------------------------------------------------------

import json
import numpy as np

# ----------------------------------------------------------
# Local imports
# ----------------------------------------------------------

# import different maps into this module
from . import mapping
from . import utils
from . import dataframe
from . import sample

# ----------------------------------------------------------
# Exceptions
# ----------------------------------------------------------

class LoadingException(Exception):
    """Error when loading Genotype Phenotype map data. """

# ----------------------------------------------------------
# Base class for constructing a genotype-phenotype map
# ----------------------------------------------------------

class GenotypePhenotypeMap(dataframe.GenotypePhenotypeDataFrame, mapping.BaseMap):
    """Main object for containing genotype-phenotype map data. Efficient memory storage,
    fast-ish mapping, graphing, and simulations.

    Parameters
    ----------
    wildtype : string
        wildtype sequence.
    genotypes : array-like
        list of all genotypes in system. Must be a complete system.
    phenotypes : array-like
        List of phenotypes in the same order as genotypes.
    mutations : dict
        Dictionary that maps each site indice to their possible substitution alphabet.
    n_replicates : int
        number of replicate measurements comprising the mean phenotypes
    include_binary : bool (default=True)
        Construct a binary representation of the space.
    """
    # ------------------------------------------------------------
    # Hidden methods for mapping object
    # ------------------------------------------------------------

    def sample(self, n_samples=1, genotypes=None, fraction=1.0, derived=True):
        """Generate artificial data sampled from phenotype and percent error.

        Parameters
        ----------
        n_samples : int
            Number of samples to take from space

        fraction : float
            fraction of space to sample.

        Returns
        -------
        samples :  Sample object
            returns this object with all stats on experiment
        """
        if genotypes is None:
            # make sure fraction is float between 0 and 1
            if fraction < 0 or fraction > 1:
                raise Exception("fraction is invalid.")
            # fractional length of space.
            frac_length = int(fraction * self.n)
            # random genotypes and phenotypes to sample
            random_indices = np.sort(np.random.choice(range(self.n), size=frac_length, replace=False))
            # If sample must include derived, set the last random_indice to self.n-1
            if derived:
                random_indices[-1] = self.n-1
        else:
            # Mapping from genotypes to indices
            mapping = self.map("genotypes", "indices")
            # Construct an array of genotype indices to sample
            random_indices = [mapping[g] for g in genotypes]

        # initialize arrays
        phenotypes = np.empty((len(random_indices), n_samples), dtype=float)
        genotypes = np.empty((len(random_indices), n_samples), dtype='<U'+str(self.length))

        # If errors are present, sample from error distribution
        try:
            # Iterate through "seen" genotypes and sample from their distributions
            for i in range(len(random_indices)):
                index = random_indices[i]
                seq = self.genotypes[index]
                # Build genotype array
                genotypes[i] = np.array([seq for j in range(n_samples)])
                stdevs = self.err.upper
                phenotypes[i] = stdevs[index] * np.random.randn(n_samples) + self.phenotypes[index]
        except AttributeError:
            # Can't sample if no error distribution is given.
            if n_samples != 1:
                raise Exception("Won't create samples if sample error is not given.")
            genotypes = np.array([self.genotypes[i] for i in random_indices])
            phenotypes = np.array([self.phenotypes[i] for i in random_indices])

        # Create a sample object
        samples = sample.Sample(self, genotypes, phenotypes, random_indices)
        return samples

    def subspace(self, genotype1, genotype2=None, mutations=None):
        """Creates a sub-GenotypePhenotypeMap object between two genotypes.

        Parameters
        ----------
        genotype1 :
            retu
        genotype2 :
            Farthest mutant sequence. Will assume from mutations dictionary if None
        mutations : dict
            Mutations dictionary

        Returns
        -------
        GenotypePhenotypeMap: obj
            Sub-GenotypePhenotypeMap between two sequences.
        """

        if genotype2 != None and len(genotype1) != len(genotype2):
            err = "genotypes must have the same length to calculate subspace.\n"
            raise ValueError(err)

        # Construct the mutations dictionary if not given (assumes binary)
        if mutations is None:
            mutations = utils.binary_mutations_map(genotype1, genotype2)
        # Construct binary encoding
        encoding = utils.encode_mutations(genotype1, mutations)
        # Construct the subspace
        wildtype = genotype1
        genotypes, binary = utils.construct_genotypes(encoding)
        # Get old genotype-phenotype mapping
        mapping = self.map("genotypes", "phenotypes")
        known, phenotypes = [], []
        for g in genotypes:
            try:
                phenotypes.append(mapping[g])
                known.append(g)
            except KeyError: pass
        # Get stdeviations
        if self.stdeviations is not None:
            mappingstd = self.map("genotypes", "stdeviations")
            stdeviations = [mappingstd[g] for g in known]
        else:
            stdeviations = None
        # Create GenotypePhenotypeMap object
        return GenotypePhenotypeMap(wildtype,
            known,
            phenotypes,
            stdeviations=stdeviations,
            log_transform=self.log_transform,
            mutations=mutations,
            n_replicates=self.n_replicates,
            logbase=self.logbase)
