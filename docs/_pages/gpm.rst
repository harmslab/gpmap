GenotypePhenotypeMap
====================

The ``GenotypePhenotypeMap`` class is main entry-point to the ``gpmap`` package.
It offers intuitive and useful methods and attributes
for analyzing genotype-phenotype data. We've also created a number of other packages that
easily interact with the ``GenotypePhenotypeMap``.

Example
-------

.. code-block:: python

    from seqspace import GenotypePhenotypeMap

    # Create list of genotypes and phenotypes
    wildtype = "AA"
    genotypes = ["AA", "AV", "AM", "VA", "VV", "VM"]
    phenotypes = [1.0, 1.1, 1.4, 1.5, 2.0, 3.0]

    # Create GenotypePhenotypeMap object
    gpm = GenotypePhenotypeMap(wildtype, genotypes, phenotypes)


Interface
---------

.. autoclass:: seqspace.gpm.GenotypePhenotypeMap
    :members:


BinaryMap
---------
All ``GenotypePhenotypeMap`` objects append a ``BinaryMap`` instance to a the ``binary``
attribute. The ``BinaryMap`` class creates a binary representation of all genotypes and maps
them to a genotype-phenotype map. Most attributes in the ``GenotypePhenotypeMap`` also exist in
the ``binary`` object, updated with the binary genotype representations.
