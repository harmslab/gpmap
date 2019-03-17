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
import gpmap.utils as utils
import gpmap.errors as errors


class GenotypePhenotypeMap(object):
    """Object for containing genotype-phenotype map data.

    Parameters
    ----------
    wildtype : string
        wildtype sequence.

    genotypes : array-like
        list of all genotypes in system. Must be a complete system.

    phenotypes : array-like
        List of phenotypes in the same order as genotypes.

    mutations : dict
        Dictionary that maps each site indice to their possible substitution
        alphabet.

    n_replicates : int
        number of replicate measurements comprising the mean phenotypes

    include_binary : bool (default=True)
        Construct a binary representation of the space.

    Attributes
    ----------
    data : pandas.DataFrame
        The core data object. Columns are 'genotypes', 'phenotypes',
        'n_replicates', 'stdeviations', and (option) 'binary'.

    complete_data : pandas.DataFrame (optional, created by BinaryMap)
        A dataframe mapping the complete set of genotypes possible, given
        the mutations dictionary. Contains all columns in `data`. Any missing
        data is reported as NaN.

    missing_data : pandas.DataFrame (optional, created by BinaryMap)
        A dataframe containing the set of missing genotypes; complte_data -
        data. Two columns: 'genotypes' and 'binary'.

    binary : BinaryMap
        object that gives you (the user) access to the binary representation
        of the map.
    """
    def __init__(self, wildtype, genotypes, phenotypes,
                 stdeviations=None,
                 mutations=None,
                 n_replicates=1,
                 **kwargs):

        # Set mutations; if not given, assume binary space.
        if mutations is not None:
            # Make sure the keys in the mutations dict are integers, not
            # strings.
            self._mutations = dict([(int(key), val)
                                   for key, val in mutations.items()])
        else:
            # Get mutations dict from genotypes.
            mutations = utils.genotypes_to_mutations(genotypes)
            self._mutations = mutations

        # Leftover kwargs become metadata that is ignored.
        self.metadata = kwargs

        # Set wildtype.
        self._wildtype = wildtype

        # Store data in DataFrame
        data = dict(
            genotypes=genotypes,
            phenotypes=phenotypes,
            n_replicates=n_replicates,
            stdeviations=stdeviations
        )
        self.data = pd.DataFrame(data)

        # Add binary representation
        self.add_binary()

        # Construct the error maps
        self._add_error()

    def _repr_html_(self):
        """Represent the GenotypePhenotypeMap as an html table."""
        return self.data.to_html()

    def map(self, attr1, attr2):
        """Dictionary that maps attr1 to attr2."""
        return dict(zip(getattr(self, attr1), getattr(self, attr2)))

    @classmethod
    def read_dataframe(cls, dataframe, wildtype, **kwargs):
        """Construct a GenotypePhenotypeMap from a dataframe."""
        # Required arguments
        df = dataframe
        self = cls(wildtype,
                   df.genotypes,
                   df.phenotypes,
                   stdeviations=df.stdeviations,
                   n_replicates=df.n_replicates,
                   **kwargs)
        return self

    @classmethod
    def read_pickle(cls, filename, **kwargs):
        """Read GenotypePhenotypeMap from pickle"""
        with open(filename, 'rb') as f:
            self = pickle.load(f)

        # Check input
        if not isinstance(self, GenotypePhenotypeMap):
            raise Exception("Pickle file does not contain a GenotypePhenotypeMap.")

        return self

    @classmethod
    def read_csv(cls, fname, wildtype, **kwargs):
        """"""
        dtypes = dict(
            genotypes=str,
            phenotypes=float,
            stdeviations=float,
            n_replicates=int
        )
        df = pd.read_csv(fname, dtype=dtypes)
        self = cls.read_dataframe(df, wildtype, **kwargs)
        return self

    @classmethod
    def read_excel(cls, fname, wildtype, **kwargs):
        """"""
        dtypes = dict(
            genotypes=str,
            phenotypes=float,
            stdeviations=float,
            n_replicates=int
        )
        df = pd.read_excel(fname, dtype=dtypes)
        self = cls.read_dataframe(df, wildtype, **kwargs)
        return self

    @classmethod
    def read_json(cls, filename, **kwargs):
        """Load a genotype-phenotype map directly from a json file.
        The JSON metadata must include the following attributes

        Note
        ----
        Keyword arguments override input that is loaded from the JSON file.
        """
        # Open, json load, and close a json file
        with open(filename, "r") as f:        
            metadata = json.load(f)
            data = metadata["data"]

        if "wildtype" in metadata:
            wildtype = metadata["wildtype"]
            metadata.pop("wildtype")

        # Check keys in dictionary.
        if not all(key in data for key in ["genotypes", "phenotypes", "stdeviations", "n_replicates"]):
            raise Exception('The "data" field must have the following keys: '
                            'genotypes", "phenotypes", "stdeviations", "n_replicates"')

        # Create an instance
        gpm = cls(
            wildtype,
            data["genotypes"],
            data["phenotypes"],
            stdeviations=data["stdeviations"],
            n_replicates=data["n_replicates"],
            **metadata
        )
        return gpm

    @classmethod
    def from_dict(self, data):
        """"""

    # ----------------------------------------------------------
    # Writing methods
    # ----------------------------------------------------------

    def to_pickle(self, filename, **kwargs):
        """Write GenotypePhenotypeMap object to a pickle file.
        """
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    def to_excel(self, filename, **kwargs):
        """Write genotype-phenotype map to excel spreadsheet.

        Keyword arguments are passed directly to Pandas dataframe to_excel
        method.

        Parameters
        ----------
        filename : str
            Name of file to write out.
        """
        self.data.to_excel(filename, **kwargs)

    def to_csv(self, filename, **kwargs):
        """Write genotype-phenotype map to csv spreadsheet.

        Keyword arguments are passed directly to Pandas dataframe to_csv
        method.

        Parameters
        ----------
        filename : str
            Name of file to write out.
        """
        self.data.to_csv(filename, **kwargs)

    def to_dict(self, complete=False):
        """Write genotype-phenotype map to dict."""
        if complete:
            data = self.complete_data.to_dict('list')
        else:
            data = self.data.to_dict('list')

        metadata = {
            "wildtype": self.wildtype,
            "mutations": self.mutations,
            "data": data
        }
        metadata.update(**self.metadata)
        return metadata

    def to_json(self, filename, complete=False):
        """Write genotype-phenotype map to json file. If no filename is given
        returns
        """
        # Get metadata.
        data = self.to_dict(complete=complete)

        # Write to json file.
        with open(filename, "w") as f:
            json.dump(data, f)

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

    @wildtype.setter
    def wildtype(self, wildtype):
        """If a wildtype is given after init, rebuild binary genotypes."""
        self._wildtype = wildtype
        self.add_binary()

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
        return self.data.genotypes.values

    @property
    def binary(self):
        """Binary representation of genotypes."""
        return self.data.binary.values

    @property
    def phenotypes(self):
        """Get the phenotypes of the system. """
        return self.data.phenotypes.values

    @property
    def stdeviations(self):
        """Get stdeviations"""
        return self.data.stdeviations.values

    @property
    def n_replicates(self):
        """Return the number of replicate measurements made of the phenotype"""
        return self.data.n_replicates.values

    @property
    def index(self):
        """Return numpy array of genotypes position. """
        return self.data.index.values

    def _add_error(self):
        """Store error maps"""
        self.std = errors.StandardDeviationMap(self)
        self.err = errors.StandardErrorMap(self)

    def add_binary(self):
        """Build a binary representation of set of genotypes.
        """
        # Construct a binary representation for each genotype
        binary = utils.genotypes_to_binary(self.wildtype,
            self.genotypes,
            self.mutations
        )

        # Add this as a column to the map.
        self.data['binary'] = binary

    def get_missing_genotypes(self):
        """Get all genotypes missing from the complete genotype-phenotype map."""
        return utils.get_missing_genotypes(
            self.genotypes,
            mutations=self.mutations
        )

    def get_all_possible_genotypes(self):
        """Get the complete set of genotypes possible. There is no particular order
        to the genotypes. Consider sorting.
        """
        # Get all genotypes.
        return mutations_to_genotypes(self.mutations, wildtype=self.wildtype)