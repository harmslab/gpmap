import numpy as np
from .gpm import GenotypePhenotypeMap

# ----------------------------------------------------------
# Sampling from Genotype-phenotype
# ----------------------------------------------------------


def sample_genotypes(gpm, fraction=1.0):
    """Return new GenotypePhenotypeMap that's a genotype sampling of given
    GenotypePhenotypeMap. Samples without replacement.

    The fraction will round to the lowest whole number of genotypes.
    """
    # Shorten gpm.data
    d = gpm.data

    # Determine size
    n = len(d)
    size = int(fraction * n)

    # Sample genotypes
    index = np.random.choice(d.index, size=size, replace=False)

    # Return GenotypePhenotypeMap
    return GenotypePhenotypeMap(gpm.wildtype,
                                d.genotypes[index],
                                d.phenotypes[index],
                                mutations=gpm.mutations,
                                stdeviations=d.stdeviations[index],
                                n_replicates=d.n_replicates[index])

def sample_phenotypes(gpm, n_replicates=1):
    """Return new GenotypePhenotypeMap that's a phenotype sampling of given
    GenotypePhenotypeMap. Samples without replacement.

    Samples from normal distributions around each phenotype.
    """
    d = gpm.data

    # Draw samples
    samples = np.random.normal(loc=d.phenotypes,
                               scale=d.stdeviations,
                               size=(n_replicates, len(d.phenotypes)))

    phenotypes = np.mean(samples, axis=0)
    stdeviations = np.std(samples, axis=0)

    # Return GenotypePhenotypeMap
    return GenotypePhenotypeMap(gpm.wildtype,
                                d.genotypes,
                                phenotypes,
                                mutations=gpm.mutations,
                                stdeviations=stdeviations,
                                n_replicates=n_replicates)
