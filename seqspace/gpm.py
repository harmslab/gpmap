# Main mapping object to be used the epistasis models in this package.
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
from seqspace.base import BaseMap
from seqspace.binary import BinaryMap
from seqspace.mutations import MutationMap
from seqspace.raw import RawMap
from seqspace.graph import GenotypePhenotypeGraph
from seqspace.errors import (StandardDeviationMap, 
                            StandardErrorMap)

from seqspace.stats import corrected_sterror

# import utils used into module.
from seqspace.utils import (hamming_distance, 
                            binary_mutations_map, 
                            farthest_genotype, 
                            encode_mutations, 
                            construct_genotypes)


# ----------------------------------------------------------
# Exceptions
# ----------------------------------------------------------

class LoadingException(Exception):
    """ Error when loading Genotype Phenotype map data. """
    
# ----------------------------------------------------------
# Base class for constructing a genotype-phenotype map
# ----------------------------------------------------------

class GenotypePhenotypeMap(BaseMap):

    def __init__(self, wildtype, genotypes, phenotypes, 
                    stdeviations=None, 
                    log_transform=False, 
                    mutations=None, 
                    n_replicates=1):
        """
            Construct a full genotype phenotype mapping object.

            Things included:
            ---------------
            1. Binary representation of all genotypes mapped to proper phenotypes
            2. Mutations mapped to their binary encoding
            3. NetworkX graph representation.

            Args:
            ----
            genotypes: array-like
                list of all genotypes in system. Must be a complete system.
            phenotypes: array-like
                List of phenotypes in the same order as genotypes.
            log_transform: boolean (default = False)
                Set to True to log tranform the phenotypes.
            mutations: dict
                Dictionary that maps each site indice to their possible substitution alphabet.

            Returns:
            -------
            GenoPhenoMap object


            GenoPhenoMap.Mutations --> includes all mutational mapping
            GenoPhenoMap.Binary --> genotypes mapped to their binary representations of the space.

        """

        # Set mutations; if not given, assume binary space.
        if mutations is not None:
            # Make sure the keys in the mutations dict are integers, not strings.
            self.mutations = dict([(int(key), val) for key, val in mutations.items()])
        else:
            mutant = farthest_genotype(wildtype, genotypes)
            mutations = binary_mutations_map(wildtype, mutant)
            self.mutations = mutations

        # Set initial properties fo GPM
        self.wildtype = wildtype
        self.genotypes = np.array(genotypes)
        self.log_transform = log_transform
        self.phenotypes = np.array(phenotypes)
        self.n_replicates = n_replicates

        # Initialize Mutational mapping
        self.Mutations = MutationMap()
        self.Mutations.mutations = mutations
        self.Mutations.wildtype = wildtype

        # Built the binary representation of the genotype-phenotype.
        # Constructs a complete sequence space and stores genotypes missing in the
        # data as an attribute, `missing_genotypes`.
        self._construct_binary()
  
  
        # Set up the error mapping
        if stdeviations is not None:
            
            self.stdeviations = np.array(stdeviations)
        
            if self.log_transform is True:
    
                # Add errors to the Raw map
                try:
                    # Add to the Raw map
                    self.Raw.stdeviations = self.stdeviations
                    self.Raw.std = StandardDeviationMap(self.phenotypes, self.stdeviations, log_transform=False)
                    self.Raw.err = StandardErrorMap(self.phenotypes, self.stdeviations, log_transform=False, n_replicates=self.n_replicates)
                    
                except AttributeError:
                    raise Exception("A RawMap must be initialized as an attribute before we can transform the errors.")
    
            # Set up all statistics for error.
            self.std = StandardDeviationMap(self.Raw.phenotypes, self.stdeviations, log_transform=self.log_transform)
            self.err = StandardErrorMap(self.Raw.phenotypes, self.stdeviations, log_transform=self.log_transform, n_replicates=self.n_replicates)
    
            # If a binary map exists
            if hasattr(self, "Binary"):    
                # Set up all statistics for error.
                self.Binary.std = StandardDeviationMap(self.Raw.phenotypes, self.stdeviations, log_transform=self.log_transform)
                self.Binary.err = StandardErrorMap(self.Raw.phenotypes, self.stdeviations, log_transform=self.log_transform, n_replicates=self.n_replicates)
        
        else:
            self.stdeviations = None
        
        # Add a networkx graph object
        self.Graph = GenotypePhenotypeGraph(self)


    # ----------------------------------------------------------
    # Class method to load from source
    # ----------------------------------------------------------

    @classmethod
    def from_json(cls, filename, **kwargs):
        """ Load a genotype-phenotype map directly from a json file. 
        
            The JSON metadata must include the following attributes
        """
        # Open, json load, and close a json file
        f = open(filename, "r")
        data = json.load(f)
        f.close()
        
        # Grab all properties from data-structure
        args = ["wildtype", "genotypes", "phenotypes"]
        options = {
            "stdeviations": None,
            "log_transform": False, 
            "mutations": None,
            "n_replicates": 1,
        }
        
        # Grab all arguments and order them
        for i in range(len(args)):
            # Get all args
            try: 
                args[i] = data[args[i]]
            except KeyError:
                raise LoadingException("""No `%s` property in json data. Must have %s for GPM construction. """ % (args[i], args[i]) )
        
        # Get all options for map and order them
        for key in options:
            # See if options are in json data
            try:
                options[key] = data[key]
            except:
                pass
        
        # Override any properties with specified kwargs passed directly into method
        options.update(kwargs)
        
        # Create an instance
        gpm = cls(args[0], args[1], args[2], **options)
        return gpm
    
    # ----------------------------------------------------------
    # Properties of the map
    # ----------------------------------------------------------


    @property
    def length(self):
        """ Get length of the genotypes. """
        return self._length

    @property
    def n(self):
        """ Get number of genotypes, i.e. size of the genotype-phenotype map. """
        return self._n

    @property
    def log_transform(self):
        """ Boolean argument telling whether space is log transformed. """
        return self._log_transform

    @property
    def wildtype(self):
        """ Get reference genotypes for interactions. """
        return self._wildtype

    @property
    def mutations(self):
        """ Get the furthest genotype from the wildtype genotype. """
        return self._mutations

    @property
    def genotypes(self):
        """ Get the genotypes of the system. """
        return self._genotypes

    @property
    def missing_genotypes(self):
        """ Genotypes that are missing from the complete genotype-to-phenotype map."""
        return self._missing_genotypes

    @property
    def complete_genotypes(self):
        """ All possible genotypes in the complete genotype space"""
        return np.concatenate((self.genotypes, self.missing_genotypes))

    @property
    def phenotypes(self):
        """ Get the phenotypes of the system. """
        return self._phenotypes
    
    @property
    def n_replicates(self):
        """Return the number of replicate measurements made of the phenotype"""
        return self._n_replicates

    @property
    def indices(self):
        """ Return numpy array of genotypes position. """
        return self._indices

    # ----------------------------------------------------------
    # Setter methods
    # ----------------------------------------------------------

    @log_transform.setter
    def log_transform(self, boolean):
        """ True/False to log transform the space. """
        self._log_transform = boolean

    @genotypes.setter
    def genotypes(self, genotypes):
        """ Set genotypes from ordered list of sequences. """
        self._n = len(genotypes)
        self._length = len(genotypes[0])
        self._genotypes = np.array(genotypes)
        self._indices = np.arange(self.n)

    @wildtype.setter
    def wildtype(self, wildtype):
        """ Set the reference genotype among the mutants in the system. """
        self._wildtype = wildtype
        self.Mutations.wildtype = wildtype

    @mutations.setter
    def mutations(self, mutations):
        """ Set the mutation alphabet for all sites in wildtype genotype.

            `mutations = { site_number : alphabet }`. If the site
            alphabet is note included, the model will assume binary
            between wildtype and derived.

            ```
            mutations = {
                0: [alphabet],
                1: [alphabet],

            }
            ```

        """
        if type(mutations) != dict:
            raise TypeError("mutations must be a dict")
        self._mutations = mutations
        self.Mutations = MutationMap()
        self.Mutations.mutations = mutations
        self.Mutations.n = len(mutations)

    @phenotypes.setter
    def phenotypes(self, phenotypes):
        """ Set phenotypes from ordered list of phenotypes

            Args:
            -----
            phenotypes: array-like or dict
                if array-like, it musted be ordered by genotype; if dict,
                this method automatically orders the phenotypes into numpy
                array.
        """
        if type(phenotypes) is dict:
            _phenotypes = self._if_dict(phenotypes)
        else:
            if len(phenotypes) != len(self._genotypes):
                raise ValueError("Number of phenotypes does not equal number of genotypes.")
            else:
                _phenotypes = phenotypes

        # log transform if log_transform = True. Raw phenotypes are stored in an separate object
        if self.log_transform is True:
            self.Raw = RawMap()
            self.Raw.phenotypes = _phenotypes
            _phenotypes = np.log10(_phenotypes)
            
        self._phenotypes = _phenotypes

        # Set binary phenotypes if binary exists... assumes
        # that binary sequences are sorted to match raw genotypes.
        if hasattr(self, "Binary"):
            self.Binary.phenotypes = phenotypes


    @n_replicates.setter
    def n_replicates(self, n_replicates):
        """Set the number of replicate measurements taken of phenotypes"""
        self._n_replicates = n_replicates


    # ------------------------------------------------------------
    # Hidden methods for mapping object
    # ------------------------------------------------------------

    def _construct_binary(self):
        """ Encode the genotypes an ordered binary set of genotypes with
            wildtype as reference state (ref is all zeros).

            This method maps each genotype to their binary representation
            relative to the 'wildtype' sequence.
        """
        # Initialize Binary representation of genotype-phenotype map
        self.Binary = BinaryMap()

        # Encode mutations as individual binary representation
        self.Binary.encoding = encode_mutations(self.wildtype, self.mutations)

        # Use encoding map to construct binary presentation for any type of alphabet
        unsorted_genotypes, unsorted_binary = construct_genotypes(self.Binary.encoding)

        # length of binary strings
        self.Binary._length = len(unsorted_binary[0])

        # Sort binary representation to match genotypes
        geno2index = self.get_map("genotypes", "indices")
        binary = np.empty(self.n, dtype="<U" + str(self.Binary.length))

        # Sort the genotypes by looking for them in the data.
        
        # Sort the missing genotypes into a separate object
        missing_genotypes = list()
        missing_binary = list()
        
        for i in range(len(unsorted_genotypes)):
            
            # Keep and sort genotype if it exists in data.
            try:
                
                index = geno2index[unsorted_genotypes[i]]
                binary[index] = unsorted_binary[i]
            
            except KeyError:
                
                missing_genotypes.append(unsorted_genotypes[i])
                missing_binary.append(unsorted_binary[i])

        # Set the missing genotypes
        self._missing_genotypes = np.array(missing_genotypes)
        self.Binary._missing_genotypes= np.array(missing_binary)

        # Set binary attributes to sorted genotypes
        self.Binary.genotypes = binary
        self.Binary.indices = self.indices[:]

        # Grab phenotypes if given by user. Otherwise, pass.
        try:
            
            self.Binary.phenotypes = self.phenotypes[:]
        
        except AttributeError:
            
            pass
