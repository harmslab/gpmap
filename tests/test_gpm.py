import pytest

from gpmap import utils
from gpmap import GenotypePhenotypeMap

import numpy as np

import time, os

def _compare_gpmap(gpm1,gpm2):
    """
    Compare key features of two genotype-phenotype maps.
    """

    array_features = ["genotypes","phenotypes","mutations","binary","index",
                      "stdeviations"]
    for a in array_features:
        assert np.array_equal(gpm1.__getattribute__(a),
                              gpm2.__getattribute__(a))

    general_features = ["wildtype","mutant","n","length"]
    for g in general_features:
        assert gpm1.__getattribute__(g) == gpm2.__getattribute__(g)

def _compare_gmap_to_input_dict(gpm,input_dict):
    """
    Compare a bunch of features of a genotype-phenotype map to data from a
    pytest input dictionary used to create it.
    """

    array_features = ["genotypes","phenotypes","mutations","binary","stdeviations"]
    for a in array_features:
        gpm_value = gpm.__getattribute__(a)

        try:
            dict_value = input_dict[a]
        except KeyError:
            continue

        assert len(gpm_value) == len(dict_value)
        for i in range(len(gpm_value)):
            assert gpm_value[i] == dict_value[i]

    general_features = ["wildtype","mutant","n","length"]
    for g in general_features:
        gpm_value = gpm.__getattribute__(g)
        try:
            dict_value = input_dict[g]
        except KeyError:
            continue

        assert gpm_value == dict_value


def test_constructor(binary_test_data):
    """
    Test basic use of GenotypePhenotypeMap constructor.
    """

    # Constructor call
    gpm = GenotypePhenotypeMap(wildtype=binary_test_data["wildtype"],
                               genotypes=binary_test_data["genotypes"],
                               log_transform=False)

    assert isinstance(gpm,GenotypePhenotypeMap)

    stripped_dict = {"wildtype":binary_test_data["wildtype"],
                     "genotypes":binary_test_data["genotypes"],
                     "binary":binary_test_data["binary"],
                     "mutant":binary_test_data["mutant"],
                     "mutations":binary_test_data["mutations"]}

    _compare_gmap_to_input_dict(gpm,stripped_dict)


def test_nonbinary_constructor(mixed_test_data):
    """
    Test GenotypePhenotypeMap constructor on non-binary input data.
    """

    # Constructor call
    gpm = GenotypePhenotypeMap(wildtype=mixed_test_data["wildtype"],
                               genotypes=mixed_test_data["genotypes"],
                               phenotypes=mixed_test_data["phenotypes"],
                               log_transform=False)

    assert isinstance(gpm,GenotypePhenotypeMap)

    _compare_gmap_to_input_dict(gpm,mixed_test_data)


def test_read_json(test_json):
    """
    Test basic reading from json.
    """
    read_gpm = GenotypePhenotypeMap.read_json(test_json)
    assert isinstance(read_gpm, GenotypePhenotypeMap)


def test_read_csv(test_csv):
    """
    Test basic reading from csv.
    """
    read_gpm = GenotypePhenotypeMap.read_csv(test_csv, wildtype='AAA')
    assert isinstance(read_gpm, GenotypePhenotypeMap)


def test_data_integrity_csv(binary_test_data,tmp_path):
    """
    Write to a csv and make sure it reads back in properly.
    """

    gpm = GenotypePhenotypeMap(wildtype=binary_test_data["wildtype"],
                               genotypes=binary_test_data["genotypes"],
                               phenotypes=binary_test_data["phenotypes"],
                               stdeviations=binary_test_data["errors"],
                               log_transform=False)

    # Write out a csv file
    out_file = os.path.join(tmp_path,"tmp.csv")
    gpm.to_csv(out_file)
    assert os.path.exists(out_file)

    # Should fail without wildtype specified
    with pytest.raises(TypeError):
        gpm_read = GenotypePhenotypeMap.read_csv(out_file)

    gpm_read = GenotypePhenotypeMap.read_csv(out_file,wildtype="AAA")

    # Make sure the written and read gpmaps ar ethe same
    _compare_gpmap(gpm,gpm_read)


#def test_std(tmp_gpm_file):
#    """Test raw errors"""
#    np.testing.assert_array_equal(tmp_gpm_file.std.upper, tmp_gpm_file.stdeviations)


def test_missing_genotypes(binary_test_data):
    """
    Test that missing genotypes are identified.
    """

    gpm = GenotypePhenotypeMap(wildtype=binary_test_data["wildtype"],
                               genotypes=binary_test_data["genotypes"],
                               phenotypes=binary_test_data["phenotypes"],
                               log_transform=False)


    # Choose a subset indices for testing missing genotypes.
    index = np.arange(0, 8)
    np.random.shuffle(index)
    # Get subset genotypes
    chosen_g = np.array(gpm.genotypes)[index[:4]]
    missing_g = np.array(gpm.genotypes)[index[4:]]

    # Get subset phenotypes
    chosen_p = np.array(gpm.phenotypes)[index[:4]]

    gpm_missing_g = GenotypePhenotypeMap(wildtype=gpm.wildtype,
                                         genotypes=chosen_g,
                                         phenotypes=chosen_p,
                                         log_transform=False,
                                         mutations=gpm.mutations)

    assert np.array_equal(gpm_missing_g.genotypes, chosen_g)
    assert np.array_equal(np.sort(gpm_missing_g.get_missing_genotypes()),
                          np.sort(missing_g))
