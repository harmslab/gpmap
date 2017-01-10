GenotypePhenotypeMap
====================

The ``GenotypePhenotypeMap`` class is the primary tool provided by the ``gpmap`` package.
It creates intuitive and useful mapping on the fly. It appends methods and attributes
to make analyzing genotype-phenotype data easy. We've create other packages that
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
Attached to the GenotypePhenotypeMap is a ``BinaryMap`` object.

The ``BinaryMap`` class creates a binary representation of all genotypes and maps
them to a genotype-phenotype map. Binary representations are useful for many reasons,
like modeling evolutionary paths and analyzing epistatic interactions (see epistasis_)

All ``GenotypePhenotypeMap`` objects append a ``BinaryMap`` instance to a the ``binary``
attribute. Most attributes in the ``GenotypePhenotypeMap`` also exist in under
the ``binary`` attribute, updated with the binary genotype representations.

.. _epistasis: http://epicstasis.readthedocs.io/
