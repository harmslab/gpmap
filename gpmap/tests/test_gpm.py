import numpy as np
import pytest
import time

from gpmap import utils
from ..gpm import GenotypePhenotypeMap


@pytest.fixture()
def tmp_gpm_file():
    wildtype = "AAA"

    genotypes = [
        "AAA",
        "AAB",
        "ABA",
        "BAA",
        "ABB",
        "BAB",
        "BBA",
        "BBB"
    ]

    binary = [
        '000',
        '001',
        '010',
        '100',
        '011',
        '101',
        '110',
        '111'
    ]

    mutations = {
        0: ["A", "B"],
        1: ["A", "B"],
        2: ["A", "B"],
    }
    phenotypes = np.random.rand(len(genotypes))
    tmp_gpm_file = GenotypePhenotypeMap(wildtype=wildtype,
                                        genotypes=genotypes,
                                        phenotypes=phenotypes,
                                        log_transform=False,
                                        mutations=mutations)
    return tmp_gpm_file


def test_read_json(tmp_gpm_file):
    """Test reading from json"""
    read_gpm = GenotypePhenotypeMap.read_json("data/test_data.json")
    # Test instance was created
    assert isinstance(read_gpm, GenotypePhenotypeMap)
    # Test elements align
    np.testing.assert_array_equal(tmp_gpm_file.genotypes, read_gpm.genotypes)


def test_read_csv():
    """Test reading from csv"""
    read_gpm = GenotypePhenotypeMap.read_csv(fname="data/test_data.csv", wildtype='AAA')
    # Test instance was created
    assert isinstance(read_gpm, GenotypePhenotypeMap)


def test_data_integrity_csv(tmp_gpm_file):
    # Export tmp_gpm_file to csv, give it time to generate file, and then check integrity
    # of the data read by "read_csv". This is necessary because phenotypes are randomly generated.
    tmp_gpm_file.to_csv(filename="data/test_data.csv")
    time.sleep(1)
    read_gpm = GenotypePhenotypeMap.read_csv(fname="data/test_data.csv", wildtype='AAA')
    # Test elements align
    np.testing.assert_array_equal(tmp_gpm_file.genotypes, read_gpm.genotypes)
    np.testing.assert_array_equal(tmp_gpm_file.phenotypes, read_gpm.phenotypes)
    np.testing.assert_array_equal(tmp_gpm_file.mutations, read_gpm.mutations)
    np.testing.assert_array_equal(tmp_gpm_file.binary, read_gpm.data.binary)


def test_length(tmp_gpm_file):
    """Test genotype length."""
    assert tmp_gpm_file.length == 3


def test_n(tmp_gpm_file):
    """Test size"""
    assert tmp_gpm_file.n == 8


def test_std(tmp_gpm_file):
    """Test raw errors"""
    np.testing.assert_array_equal(tmp_gpm_file.std.upper, tmp_gpm_file.stdeviations)


def test_mutant(tmp_gpm_file):
    """Test mutant"""
    assert tmp_gpm_file.mutant == "BBB"


def test_missing_genotypes(tmp_gpm_file):
    """Test that missing genotypes are identified."""
    # Choose a subset indices for testing missing genotypes.
    index = np.arange(0, 8)
    np.random.shuffle(index)
    # Get subset genotypes
    chosen_g = np.array(tmp_gpm_file.genotypes)[index[:4]]
    missing_g = np.array(tmp_gpm_file.genotypes)[index[4:]]

    # Get subset phenotypes
    chosen_p = np.array(tmp_gpm_file.phenotypes)[index[:4]]

    gpm_missing_g = GenotypePhenotypeMap(wildtype=tmp_gpm_file.wildtype,
                                         genotypes=chosen_g,
                                         phenotypes=chosen_p,
                                         log_transform=False,
                                         mutations=tmp_gpm_file.mutations
                                         )

    np.testing.assert_array_equal(gpm_missing_g.genotypes, chosen_g)
    np.testing.assert_array_equal(np.sort(gpm_missing_g.get_missing_genotypes()), np.sort(missing_g))
