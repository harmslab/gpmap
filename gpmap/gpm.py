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
            self.add_binary(self.wildtype)

        # Construct the error maps
        self._add_error()

    @classmethod
    def read_dataframe(cls, wildtype, dataframe, **kwargs):
        """Construct a GenotypePhenotypeMap from a dataframe."""
        # Required arguments
        df = DataFrames
        genotypes = df["genotypes"]
        phenotypes = df["phenotypes"]

        # Search for optional columns
        other_items = {}
        for key in OPTIONAL_KEYS:
            try:
                other_items[key] = df[key]
            except KeyError:
                pass

        # Initialize object.
        self = cls(wildtype, genotypes, phenotypes, **other_items)
        return self

    @classmethod
    def read_excel(cls, wildtype, fname, **kwargs):
        """"""
        df = pd.read_excel(fname)
        self = cls.read_dataframe(wildtype, df)
        return self

    @classmethod
    def read_csv(cls, fname, **kwargs):
        """"""
        df = pd.read_csv(fname)
        self = cls.read_dataframe(wildtype, df)
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
    def mutations(self):
        """Get the furthest genotype from the wildtype genotype."""
        return self._mutations

    @property
    def genotypes(self):
        """Get the genotypes of the system."""
        return self.data.genotypes.values

    @property
    def missing_genotypes(self):
        """Genotypes that are missing from the complete genotype-to-phenotype
        map."""
        return self.missing_data.genotypes.values

    @property
    def complete_genotypes(self):
        """Array of sorted genotypes for the complete genotype space encoded by
        the mutations dictionary.

        **NOTE** Can only be set by the BinaryMap object.
        """
        try:
            return self.complete_data.genotypes.values
        except AttributeError:
            raise AttributeError("Looks like a BinaryMap has not been built "
                                 "yet for this map. Do this before asking for "
                                 "the complete_genotypes.")

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

    def add_binary(self, wildtype):
        """Add a BinaryMap to the GenotypePhenotypeMap. The wildtype determines
        the encoding pattern. Wildtype sites are represented as 0's.
        """
        self.binary = binary.BinaryMap(self, wildtype)
