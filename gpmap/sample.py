import numpy as np
from . import gpm
# ----------------------------------------------------------
# Sampling from Genotype-phenotype
# ----------------------------------------------------------

class Sample(object):

    def __init__(self, gpm, replicate_genotypes, replicate_phenotypes, indices=None):
        """Sample from simulated experiment """
        self._gpm = gpm
        self.replicate_genotypes = replicate_genotypes
        self.replicate_phenotypes = replicate_phenotypes
        self.n_replicates = replicate_genotypes.shape[1]
        try:
            self.genotypes = self.replicate_genotypes[:,0]
            self.phenotypes = np.mean(self.replicate_phenotypes, axis=1)
            self.stdeviations = np.std(self.replicate_phenotypes, ddof=1, axis=1)
        except IndexError:
            self.genotypes = self.replicate_genotypes
            self.phenotypes = self.replicate_phenotypes
            self.stdeviations = self.replicate_phenotypes
        self.indices = indices

    def get_gpm(self):
        """Return a Genotype-phenotype object from sample. """
        return gpm.GenotypePhenotypeMap(self._gpm.wildtype, self.genotypes, self.phenotypes,
                stdeviations=self.stdeviations,
                mutations=self._gpm.mutations,
                n_replicates=self.n_replicates)
