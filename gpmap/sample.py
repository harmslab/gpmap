import numpy as np
from .gpm import GenotypePhenotypeMap

# ----------------------------------------------------------
# Sampling from Genotype-phenotype
# ----------------------------------------------------------


def sample_genotypes(gpm, fraction=1.0, size=None):
    """Return new GenotypePhenotypeMap that's a genotype sampling of given
    GenotypePhenotypeMap. Samples without replacement.

    The fraction will round to the lowest whole number of genotypes.
    If size is given, fraction is ignored
    """
    # Shorten gpm.data
    d = gpm.data

    # Determine size
    n = len(d)
    if size is None:
        size = int(fraction * n)
    elif size > n:
        Exception('Size exceeds the number of genotypes.')

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


def subset(gpm, **rows_to_select):
    """Create a subset GenotypePhenotypeMap from a given GenotypePhenotypeMap
    using an attribute of the Map object and a list of values.

    For example, we can subset a map using the genotypes attribute and give
    a list of genotypes to include in the subset

    .. code-block:: python

        gpm2 = subset(gpm, genotypes=['AA', 'AV'])


    Alternatively, we can subset the map using the binary attribute the same
    way.

    .. code-block:: python

        gpm2 = subset(gpm, binary=['00', '01'])
    """
    d = gpm.data
    items = list(rows_to_select.items())

    # Check kwargs is proper
    if len(items) > 1:
        raise Exception('This function can only call one column.')
    elif items[0][0] not in ['genotypes', 'binary', 'index']:
        raise Exception('Attribute given not in genotype-phenotype map.')
    elif items[0][0] == 'index':
        column = items[0][1]
        d_ = d.loc[items[0][1]]
    else:
        label = items[0][0]
        column = items[0][1]
        d_ = d.loc[d[label].isin(column)]

    # Not all the
    if len(d_) != len(column):
        raise Exception("Not all rows were found.")

    # Return GenotypePhenotypeMap
    return GenotypePhenotypeMap(gpm.wildtype,
                                d_.genotypes,
                                d_.phenotypes,
                                mutations=gpm.mutations,
                                stdeviations=d_.stdeviations,
                                n_replicates=d_.n_replicates)
