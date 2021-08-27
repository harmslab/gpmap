__description__ = \
"""
Configure test environment for gpmap.
"""

import pytest

import os

@pytest.fixture(scope="module")
def test_csv():
    """
    Example csv file with genotype-phenotype map
    """

    # Find directory with test files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.abspath(os.path.join(base_dir,"data"))

    return os.path.join(data_dir,"test_data.csv")

@pytest.fixture(scope="module")
def test_json():
    """
    Example json file with genotype-phenotype map
    """

    # Find directory with test files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.abspath(os.path.join(base_dir,"data"))

    return os.path.join(data_dir,"test_data.json")

@pytest.fixture(scope="module")
def test_csv_five():
    """
    Example csv file with 5 binary sites genotype-phenotype map
    """

    # Find directory with test files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.abspath(os.path.join(base_dir,"data"))

    return os.path.join(data_dir,"test_data_five.csv")

@pytest.fixture(scope="module")
def binary_test_data():
    """
    Simple binary map.
    """

    wildtype = "AAA"

    genotypes = ["AAA", "AAB", "ABA", "ABB",
                 "BAA", "BAB", "BBA", "BBB"]

    binary = ["000", "001", "010", "011",
              "100", "101", "110", "111"]

    mutations = {0:["A", "B"],
                 1:["A", "B"],
                 2:["A", "B"]}

    phenotypes = [0.2611278, 0.60470609, 0.13114308, 0.76428437,
                  0.5018751, 0.18654072, 0.88086482, 0.18263346]

    errors = [0.05, 0.05, 0.05, 0.05,
              0.05, 0.05, 0.05, 0.05]

    return {"genotypes":genotypes,
            "binary":binary,
            "phenotypes":phenotypes,
            "errors":errors,
            "wildtype":wildtype,
            "mutations":mutations,
            "length":8,
            "n":8,
            "mutant":"BBB"}

@pytest.fixture(scope="module")
def mixed_test_data():
    """
    Map that is binary at two sites and trinary at a third site.
    """

    wildtype = "AAA"

    genotypes = ["AAA", "AAB", "AAC", "ABA",
                 "ABB", "ABC", "BAA", "BAB",
                 "BAC", "BBA", "BBB", "BBC"]

    binary = ['0000', '0010', '0001', '0100',
              '0110', '0101', '1000', '1010',
              '1001', '1100', '1110', '1101']

    mutations = {0:["A", "B"],
                 1:["A", "B"],
                 2:["A", "B","C"]}

    phenotypes = [0.60371285, 0.10893567, 0.49704416, 0.34674266,
                  0.26102007, 0.02631915, 0.44587924, 0.31596652,
                  0.87037953, 0.95649285, 0.39668621, 0.66987709]

    errors = [0.05, 0.05, 0.05, 0.05,
              0.05, 0.05, 0.05, 0.05,
              0.05, 0.05, 0.05, 0.05]

    return {"genotypes":genotypes,
            "binary":binary,
            "phenotypes":phenotypes,
            "errors":errors,
            "wildtype":wildtype,
            "mutations":mutations}
