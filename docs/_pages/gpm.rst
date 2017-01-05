GenotypePhenotypeMap
====================

The ``GenotypePhenotypeMap`` class is the primary tool provided by the ``seqspace`` package.
It creates intuitive and useful mapping on the fly. It appends methods and attributes
for analysis.

BinaryMap
---------
Attached to the GenotypePhenotypeMap is a ``BinaryMap`` object.

The ``BinaryMap`` class creates a binary representation of all genotypes and maps
them to a genotype-phenotype map. Binary representations are useful for many reasons,
like modeling evolutionary paths and analyzing epistatic interactions (see epistasis_)

All ``GenotypePhenotypeMap`` objects append a ``BinaryMap`` instance to a the ``binary``
attribute. Most attributes in the ``GenotypePhenotypeMap`` also exist in under
the ``binary`` attribute, updated with the binary genotype representations.

.. _epistasis: http://epistasis.readthedocs.io/
