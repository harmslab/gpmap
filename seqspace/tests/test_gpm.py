import numpy as np
import nose.tools as tools


from . import base
from ..gpm import GenotypePhenotypeMap

class testGenotypePhenotypeMap(base.BaseTestClass):

    def _test_attributes(self, gpm):
        # Test instance was created
        tools.assert_is_instance(gpm, GenotypePhenotypeMap)
        # Test elements align
        np.testing.assert_array_equal(gpm.genotypes, self.genotypes)
        np.testing.assert_array_equal(gpm.phenotypes, np.log10(self.phenotypes))

    def test_init(self):
        """Test initialization of class"""
        gpm = GenotypePhenotypeMap(
            self.wildtype,
            self.genotypes,
            self.phenotypes,
            log_transform=self.log_transform,
            n_replicates=self.n_replicates,
            mutations=self.mutations
        )
        # Test instance was created
        tools.assert_is_instance(gpm, GenotypePhenotypeMap)
        # Test elements align
        np.testing.assert_array_equal(gpm.genotypes, self.genotypes)
        np.testing.assert_array_equal(gpm.phenotypes, np.log10(self.phenotypes))

    def test_from_json(self):
        """Test reading from json"""
        gpm = GenotypePhenotypeMap.from_json("data.json")
        # Test instance was created
        tools.assert_is_instance(gpm, GenotypePhenotypeMap)
        # Test elements align
        np.testing.assert_array_equal(gpm.genotypes, self.genotypes)
        np.testing.assert_array_equal(gpm.phenotypes, np.log10(self.phenotypes))

    def test_length(self):
        """Test genotype length."""
        gpm = GenotypePhenotypeMap.from_json("data.json")
        tools.assert_equal(self.GPM.length, 4)

    def test_n(self):
        """Test size"""
        tools.assert_equal(self.GPM.n, 2**4)

    def test_raw_phenotypes(self):
        """Test log_transform to raw phenotypes"""
        np.testing.assert_array_equal(self.GPM.Raw.phenotypes, self.phenotypes)

    def test_missing_genotypes(self):
        """Test that missing genotypes are identified."""
        # Choose a subset indices for testing missing genotypes.
        index = np.arange(0,16)
        np.random.shuffle(index)
        # Get subset genotypes
        chosen_g = np.array(self.genotypes)[index[:8]]
        missing_g = np.array(self.genotypes)[index[8:]]

        # Get subset phenotypes
        chosen_p = np.array(self.phenotypes)[index[:8]]
        # Init GPM
        gpm = GenotypePhenotypeMap(
            self.wildtype,
            chosen_g,
            chosen_p,
            log_transform=self.log_transform,
            n_replicates=self.n_replicates,
            mutations=self.mutations,
        )
        # Test missing genotypes are returned
        np.testing.assert_array_equal(gpm.genotypes, chosen_g)
        np.testing.assert_array_equal(np.sort(gpm.missing_genotypes), np.sort(missing_g))
