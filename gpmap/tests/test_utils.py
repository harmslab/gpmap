# Import utils model.
from .. import utils

WILDTYPE = "AAA"

GENOTYPES = [
    "AAA",
    "AAB",
    "ABA",
    "BAA",
    "ABB",
    "BAB",
    "BBA",
    "BBB"
]

BINARY = [
    '000',
    '001',
    '010',
    '100',
    '011',
    '101',
    '110',
    '111'
]

MUTATIONS = {
    0: ["A", "B"],
    1: ["A", "B"],
    2: ["A", "B"],
}

def lists_are_same(list1, list2):
    """Return true if lists contain same values
    (order doesn't matter).
    """
    set1, set2 = set(list1), set(list2)

    # Find differences
    if set1.issubset(set2) and set1.issuperset(set2):
        return True

    else:
        return False


def dicts_are_same(dict1, dict2):
    """Return True if dictionaries are the same.
    """
    # Check that keys are the same.
    keys1, keys2 = list(dict1.keys()), list(dict2.keys())
    if not lists_are_same(keys1, keys2):
        return False

    # Check that values are the same.
    for key, val in dict1.items():
        if dict2[key] != val:
            return False

    return True

def test_hamming_distance():
    """Test hamming distance function."""
    s1 = "THIS IS A TEST"
    s2 = "HWIS IT A TEBT"

    assert utils.hamming_distance(s1, s2) == 4


def test_find_differences():
    """Test find differences function."""
    s1 = "THIS IS A TEST"
    s2 = "HWIS IT A TEBT"

    assert utils.find_differences(s1, s2) == [0,1,6,12]


def test_farthest_genotype():
    """Test farthest genotypes function."""
    assert utils.farthest_genotype("AAA", GENOTYPES) == "BBB"


def test_list_binary():
    """test list binary function."""
    bin = utils.list_binary(3)

    assert lists_are_same(bin, BINARY)


def test_mutations_to_encoding():
    """Test mutations to encoding."""
    encoding = utils.mutations_to_encoding(WILDTYPE, MUTATIONS)

    assert encoding[0]['A'] == '0'
    assert encoding[1]['B'] == '1'
    assert encoding[2]['A'] == '0'


def test_mutations_to_genotypes():
    """Test mutations to genotypes function."""
    genotypes = utils.mutations_to_genotypes(MUTATIONS, wildtype=WILDTYPE)

    assert lists_are_same(genotypes, GENOTYPES)


def test_genotypes_to_mutations():
    """Test genotypes to mutations function."""
    mutations = utils.genotypes_to_mutations(GENOTYPES)

    assert dicts_are_same(mutations, MUTATIONS)


def test_genotypes_to_binary():
    """Test get_missing_genotypes function."""
    encoding_table = utils.get_encoding_table(WILDTYPE, MUTATIONS)
    binary = utils.genotypes_to_binary(GENOTYPES, encoding_table)

    assert lists_are_same(binary, BINARY)


def test_get_missing_genotypes():
    """Test get_missing_genotypes function."""
    known_, missing_ = GENOTYPES[0:4], GENOTYPES[4:]

    missing = utils.get_missing_genotypes(known_, MUTATIONS)

    assert lists_are_same(missing, missing_)
