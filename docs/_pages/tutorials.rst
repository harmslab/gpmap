Tutorials
=========

Below are a set of simple tutorials for using the `seqspace` package.

GenotypePhenotypeMap
--------------------

To begin using this package, import the base class of the package, `GenotypePhenotypeMap`.

.. code-block:: python

    from gpmap import GenotypePhenotypeMap

There a few ways you initialize a GenotypePhenotypeMap object. First, you can pass
the data into the object directly:

.. code-block:: python

    wildtype = "00"
    genotypes = ["00", "01", "10", "11"]
    phenotypes = [0, 0.5, 0.5, 1]

    gpm = GenotypePhenotypeMap(wildtype, genotypes, phenotypes)


Alternatively, you can load data from a json file (using the format defined `here`_)

.. code-block:: python

    path = "data.json"
    gpm = GenotypePhenotypeMap.from_json(path)


.. _here: io.rst


Using NetworkX
--------------

.. code-block:: python

    Graph = gpm.add_networkx()

.. code-block:: python
