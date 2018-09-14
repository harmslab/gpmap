Helpful functions
=================

GPMap comes with many helpful functions for enumerating genotype-phenotype maps.
This page provides a simple list of those functions.

- :ref:`get-all-genotypes-from-mutations`
- :ref:`genotypes-to-mutations`
- :ref:`genotypes-to-binary`
- :ref:`missing-genotypes`

.. _get-all-genotypes-from-mutations:

Get all genotypes from mutations
--------------------------------

.. ipython::

  In [1]: from gpmap.utils import genotypes_to_mutations

  In [2]: genotypes = [
     ...:     "AAA",
     ...:     "AAB",
     ...:     "ABA",
     ...:     "BAA",
     ...:     "ABB",
     ...:     "BAB",
     ...:     "BBA",
     ...:     "BBB"
     ...: ]

  In [3]: genotypes_to_mutations(genotypes)
  Out[3]: {0: ['A', 'B'], 1: ['A', 'B'], 2: ['A', 'B']}



.. _`genotypes-to-mutations`:

Get mutations from a list of genotypes
--------------------------------------

.. ipython::

  In [1]: from gpmap.utils import mutations_to_genotypes

  In [2]: mutations = {0: ['A', 'B'], 1: ['A', 'B'], 2: ['A', 'B']}

  In [3]: mutations_to_genotypes(mutations)
  Out[3]: ['AAA', 'AAB', 'ABA', 'ABB', 'BAA', 'BAB', 'BBA', 'BBB']


.. _`genotypes-to-binary`:

Get binary representation of genotypes
--------------------------------------

.. ipython::

  In [1]: from gpmap.utils import genotypes_to_binary

  In [2]: wildtype = 'AAA'

  In [3]: genotypes = [
     ...:     "AAA",
     ...:     "AAB",
     ...:     "ABA",
     ...:     "BAA",
     ...:     "ABB",
     ...:     "BAB",
     ...:     "BBA",
     ...:     "BBB"
     ...: ]

  In [4]: mutations = {0: ['A', 'B'], 1: ['A', 'B'], 2: ['A', 'B']}

  In [5]: genotypes_to_binary(wildtype, genotypes, mutations)
  Out[5]: ['000', '001', '010', '100', '011', '101', '110', '111']


.. _`missing-genotypes`:

Get a list of missing genotypes from a list of genotypes
--------------------------------------------------------

.. ipython::

  In [1]: from gpmap.utils import get_missing_genotypes

  In [2]: genotypes = ["AAA","BBB"]

  In [3]: get_missing_genotypes(genotypes)
  Out[3]: ['BBA', 'BAB', 'ABB', 'ABA', 'AAB', 'BAA']
