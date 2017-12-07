import numpy as np
import pandas as pd
from .. import GenotypePhenotypeMap


def mask(gpm, mask_fraction):
    """Create a new GenotypePhenotypeMap object from a random subset of another
    GenotypePhenotypeMap.


    Returns
    -------
    true_mask_fraction : float
        the actual fraction used, since the space is discrete and likely
        won't be the exact fraction given.
    GenotypePhenotypeMap :
        the new genotype-phenotype map.
    """
    if mask_fraction > 1 or mask_fraction < 0:
        raise Exception("mask_fraction must between between 0 and 1.")

    # Calculate the number of genotypes to select
    number_to_choose = int((1 - mask_fraction) * gpm.n)

    # Calculate the true fraction (since this is a discrete space.)
    true_mask_fraction = 1 - float(number_to_choose) / gpm.n

    # Randomly choose genotypes
    index = np.random.choice(gpm.index, number_to_choose, replace=False)

    # Check n_replicates datatype
    if type(gpm.n_replicates) == int:
        n_replicates = gpm.n_replicates
    elif type(gpm.n_replicates) == pd.Series or
    type(gpm.n_replicates) == np.ndarray:
        n_replicates = gpm.n_replicates[index].reset_index(drop=True)
    else:
        raise Exception("n_replicates are not a valid dtype.")

    # return Subset genotype
    return true_mask_fraction, GenotypePhenotypeMap(
        gpm.wildtype,
        gpm.genotypes[index].reset_index(drop=True),
        gpm.phenotypes[index].reset_index(drop=True),
        mutations=gpm.mutations,
        stdeviations=gpm.stdeviations[index].reset_index(drop=True),
        n_replicates=n_replicates)
