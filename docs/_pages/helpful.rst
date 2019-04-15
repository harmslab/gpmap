Helpful functions
=================

GPMap comes with many helpful functions for enumerating genotype-phenotype maps.
This page provides a simple list of those functions.

- :ref:`get-all-genotypes-from-mutations`
- :ref:`missing-genotypes`
- :ref:`genotypes-to-mutations`
- :ref:`genotypes-to-binary`
- :ref:`get-encoding-table`

.. _get-all-genotypes-from-mutations:

Get all genotypes from mutations
--------------------------------

.. code-block:: python

  from gpmap.utils import genotypes_to_mutations

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

  mutations = genotypes_to_mutations(genotypes)  

.. _`get-encoding-table`:


Get mutation encoding table
---------------------------

.. code-block:: python

  from gpmap.utils import get_encoding_table

  wildtype = "AA"
  mutations = {
      0: ["A", "B"],
      1: ["A", "B"]
  }
  get_encoding_table(wildtype, mutations)

.. raw:: html

  <table border="1">
    <thead>
      <tr style="text-align: right;">
        <th></th>
        <th>binary_index_start</th>
        <th>binary_index_stop</th>
        <th>binary_repr</th>
        <th>genotype_index</th>
        <th>mutation_index</th>
        <th>mutation_letter</th>
        <th>wildtype_letter</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <th>0</th>
        <td>0</td>
        <td>1</td>
        <td>0</td>
        <td>0</td>
        <td>NaN</td>
        <td>A</td>
        <td>A</td>
      </tr>
      <tr>
        <th>1</th>
        <td>0</td>
        <td>1</td>
        <td>1</td>
        <td>0</td>
        <td>1</td>
        <td>B</td>
        <td>A</td>
      </tr>
      <tr>
        <th>2</th>
        <td>1</td>
        <td>2</td>
        <td>0</td>
        <td>1</td>
        <td>NaN</td>
        <td>A</td>
        <td>A</td>
      </tr>
      <tr>
        <th>3</th>
        <td>1</td>
        <td>2</td>
        <td>1</td>
        <td>1</td>
        <td>2</td>
        <td>B</td>
        <td>A</td>
      </tr>
    </tbody>
  </table>

.. _`genotypes-to-mutations`:

Get mutations from a list of genotypes
--------------------------------------

.. code-block:: python

  from gpmap.utils import mutations_to_genotypes

  mutations = {0: ['A', 'B'], 1: ['A', 'B'], 2: ['A', 'B']}

  mutations_to_genotypes(mutations)
  # ['AAA', 'AAB', 'ABA', 'ABB', 'BAA', 'BAB', 'BBA', 'BBB']


.. _`genotypes-to-binary`:

Get binary representation of genotypes
--------------------------------------

.. code-block:: python

  from gpmap.utils import genotypes_to_binary, get_encoding_table

  wildtype = 'AAA'

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

  mutations = {0: ['A', 'B'], 1: ['A', 'B'], 2: ['A', 'B']}
  table = get_encoding_table(wildtype, mutations)
  binary = genotypes_to_binary(genotypes, table)
  # ['000', '001', '010', '100', '011', '101', '110', '111']

.. _`missing-genotypes`:

Get a list of missing genotypes from a list of genotypes
--------------------------------------------------------

.. code-block:: python

  from gpmap.utils import get_missing_genotypes

  genotypes = ["AAA","BBB"]

  get_missing_genotypes(genotypes)
  # ['BBA', 'BAB', 'ABB', 'ABA', 'AAB', 'BAA']
