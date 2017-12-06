# Base class for all maps in this a genotype-phenotype map.
#
# Author: Zach Sailer
#
# -------------------------------------
# Outside imports
# -------------------------------------
import numpy as np
from collections import OrderedDict

# -------------------------------------
# Main class for building epistasis map
# -------------------------------------


class BaseMap(object):
    """Base class for all maps in this file.
    """

    def _if_dict(self, dictionary):
        """ If setter method is passed a dictionary with genotypes as keys,
        use those keys to populate array of elements in order
        """
        elements = np.empty(self._n, dtype=float)
        for i in range(self._n):
            elements[i] = dictionary[self._genotypes[i]]
        return elements

    def map(self, attr1, attr2):
        """ Return a mapping dictionary between two attributes in map.

        Parameters
        ----------
        attr1 : str
            __name__ of attribute that will be keys of dictionary
        attr2 : str
            __name__ of attribute that will the values of dictionary

        Returns
        -------
        mapping : dict
            mapping attr1 to att2
        """
        # for handling nested object attributes
        def nested_attr(main_obj, attr):
            """ Get attributes in nested objects.

            Returns the children object that is the parent object
            of the attribute of interest.

            Also returns the attribute as a string
            """
            levels = attr.split(".")
            subclasses = levels[:-1]
            attr = levels[-1]

            obj = main_obj
            for sub_obj in subclasses:
                obj = obj.__getattribute__(sub_obj)

            return obj, attr

        # Get keys
        obj, attr = nested_attr(self, attr1)
        keys = getattr(obj, attr)

        obj, attr = nested_attr(self, attr2)
        values = getattr(obj, attr)

        # Construct map
        mapping = dict(zip(keys, values))

        return mapping
