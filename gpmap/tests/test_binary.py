import numpy as np
import nose.tools as tools


from . import base
from ..gpm import GenotypePhenotypeMap
from ..errors import StandardErrorMap, StandardDeviationMap
from ..binary import BinaryMap

class testBinaryMap(base.BaseTestClass):

    def setUp(self):
        super(testBinaryMap,self).setUp()
        self.binary = self.GPM.binary

    def test_init(self):
        # Test elements align
        tools.assert_is_instance(self.binary._GPM, GenotypePhenotypeMap)
        tools.assert_is_instance(self.binary.std, StandardDeviationMap)
        tools.assert_is_instance(self.binary.err, StandardErrorMap)
        np.testing.assert_array_equal(self.binary.phenotypes, self.phenotypes)

    def test_length(self):
        tools.assert_equal(self.binary.length, 4)

    def test_phenotypes(self):
        """Test non-log_transform to raw phenotypes"""
        np.testing.assert_array_equal(self.binary.phenotypes, self.phenotypes)

    def test_std(self):
        """Test raw errors"""
        np.testing.assert_array_equal(self.binary.std.upper, self.stdeviations)
