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


class GenotypePhenotypeMap(mapping.BaseMap):
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
        the mutations dictionary. Two columns: 'genotypes' and 'binary'.

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
            phenotypes=phenotypes,
            n_replicates=n_replicates,
            stdeviations=stdeviations
        )
        self.data = pd.DataFrame(data)

        # Built the binary representation of the genotype-phenotype.
        # Constructs a complete sequence space and stores genotypes missing
        # in the data as an attribute, `missing_genotypes`.
        if include_binary:
            self.add_binary()

        # Construct the error maps
        self._add_error()

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
    def read_excel(cls, fname, wildtype, **kwargs):
        """"""
        df = pd.read_excel(fname)
        self = cls.read_dataframe(df, wildtype)
        return self

    @classmethod
    def read_csv(cls, fname, wildtype, **kwargs):
        """"""
        df = pd.read_csv(fname)
        self = cls.read_dataframe(df, wildtype)
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
        f = open(filename, "r")
        data = json.load(f)
        f.close()

        # Grab all properties from data-structure
        necessary_args = ["wildtype", "genotypes", "phenotypes"]
        options = {
            "genotypes": [],
            "phenotypes": [],
            "wildtype": [],
            "stdeviations": None,
            "mutations": None,
            "n_replicates": 1,
        }
        # Get all options for map and order them
        for key in options:
            # See if options are in json data
            try:
                options[key] = data[key]
            except KeyError:
                pass
        # Override any properties with manually entered kwargs passed directly
        # into method
        options.update(kwargs)
        args = []
        for arg in necessary_args:
            val = options.pop(arg)
            args.append(val)
        # Create an instance
        gpm = cls(args[0], args[1], args[2], **options)
        return gpm

    # ----------------------------------------------------------
    # Writing methods
    # ----------------------------------------------------------

    def to_excel(self, filename, **kwargs):
        """Write genotype-phenotype map to excel spreadsheet.

        Keyword arguments are passed directly to Pandas dataframe to_excel
        method.

        Parameters
        ----------
        filename : str
            Name of file to write out.
        """
        self.df.to_excel(filename, **kwargs)

    def to_csv(self, filename, **kwargs):
        """Write genotype-phenotype map to csv spreadsheet.

        Keyword arguments are passed directly to Pandas dataframe to_csv
        method.

        Parameters
        ----------
        filename : str
            Name of file to write out.
        """
        self.df.to_csv(filename, **kwargs)

    def to_json(self, filename):
        """Write genotype-phenotype map to json file.
        """
        # Get metadata.
        data = dict(wildtype=self.wildtype,
                    genotypes=list(self.genotypes),
                    phenotypes=list(self.phenotypes),
                    stdeviations=list(self.stdeviations),
                    mutations=self.mutations,
                    n_replicates=list(self.n_replicates.astype(float)))

        # Write to file
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
    def complete_data(self):
        """"""
        return pd.merge(self._complete_data, self.data,
                        on=['genotypes', 'binary'],
                        how='outer')
        #return pd.merge((self._complete_data, self.data),)

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
    def missing_genotypes(self):
        """Genotypes that are missing from the complete genotype-to-phenotype
        map."""
        return self.missing_data.genotypes.values

    @property
    def missing_binary(self):
        """Binary representation of missing genotypes"""
        return self.missing_data.binary.values

    @property
    def complete_genotypes(self):
        """Both missing and observed genotypes.
        """
        return self.complete_data.genotypes.values

    @property
    def complete_binary(self):
        """Complete set of genotypes as binary representation."""
        return self.complete_data.binary.values

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
        """Build a binary representation of set of
        genotypes. Also the full set of genotypes, and creates two new
        DataFrames: 'missing_data' and 'complete_data'.

        Also reindexes `data` based on `complete_data`.
        """
        # Encode mutations as binary rep.
        encoding = utils.encode_mutations(self.wildtype, self.mutations)

        # Use encoding map to construct binary presentation
        unsorted_genotypes, unsorted_binary = utils.construct_genotypes(
            encoding)

        # New data object.
        data = {'genotypes': unsorted_genotypes,
                'binary': unsorted_binary}

        # Store the complete genotype space in a DataFrame.
        self._complete_data = pd.DataFrame(data)
        self._complete_data.sort_values('genotypes', inplace=True)
        self._complete_data.reset_index(drop=True, inplace=True)

        # Mapping genotypes to index
        mapping = dict(zip(self._complete_data.genotypes,
                           self._complete_data.index))

        # Get index of observed genotypes and missing genotypes.
        observed_index = [mapping[g] for g in self.data.genotypes]
        missing_index = set(self._complete_data.index).difference(
            observed_index)

        # Reset index of main data.
        self.data.index = observed_index
        # Add a column for binary representation of genotypes.
        self.data['binary'] = pd.Series(self._complete_data.binary,
                                        index=observed_index)

        # Create a dataframe for the missing data.
        self.missing_data = pd.DataFrame(self.complete_data,
                                         index=missing_index)
