import numpy as np
import pytest

from . import base
from ..gpm import GenotypePhenotypeMap


@pytest.fixture()
def tmp_gpm_file():
    pass


class testGenotypePhenotypeMap(base.BaseTestClass):

    def _test_attributes(self, gpm):
        # Test instance was created
        assert isinstance(gpm, GenotypePhenotypeMap)
        # Test elements align
        np.testing.assert_array_equal(gpm.genotypes, self.genotypes)

    def test_init(self):
        """Test initialization of class"""
        gpm = GenotypePhenotypeMap(
            self.wildtype,
            self.genotypes,
            self.phenotypes,
            n_replicates=self.n_replicates,
            mutations=self.mutations
        )
        # Test instance was created
        tools.assert_is_instance(gpm, GenotypePhenotypeMap)
        # Test elements align
        np.testing.assert_array_equal(gpm.genotypes, self.genotypes)
        np.testing.assert_array_equal(gpm.phenotypes, self.phenotypes)

    def test_read_json(self):
        """Test reading from json"""
        gpm = GenotypePhenotypeMap.read_json("data.json")
        # Test instance was created
        tools.assert_is_instance(gpm, GenotypePhenotypeMap)
        # Test elements align
        np.testing.assert_array_equal(gpm.genotypes, self.genotypes)

    def test_length(self):
        """Test genotype length."""
        gpm = GenotypePhenotypeMap.read_json("data.json")
        tools.assert_equal(self.GPM.length, 4)

    def test_n(self):
        """Test size"""
        tools.assert_equal(self.GPM.n, 2**4)

    def test_std(self):
        """Test raw errors"""
        np.testing.assert_array_equal(self.GPM.std.upper, self.stdeviations)

    def test_mutant(self):
        """Test mutant"""
        tools.assert_equal(self.GPM.mutant, "1111")

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
