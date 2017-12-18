.. gpmap documentation master file, created by
   sphinx-quickstart on Fri Jul  8 10:41:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``gpmap``
=========

The Pandas DataFrame for genotype-phenotype (GP) map data.

.. image:: _img/gpm.png
    :align: center


The ``GenotypePhenotypeMap`` is a core object for a suite of packages written
in the `Harms Lab`_. It organizes and standardizes genotype-phenotype map data.

.. _`Harms Lab`: https://github.com/harmslab

Basic Example
-------------

.. code-block:: python

  # Import the GenotypePhenotypeMap
  from gpmap import GenotypePhenotypeMap

  # The data
  wildtype = 'AA'
  genotypes = ['AA', 'AT', 'TA', 'TT']
  phenotypes = [0.1, 0.5, 0.2, 0.8]
  stdeviations = [0.05, 0.05, 0.05, 0.05]

  # Initialize a GenotypePhenotype object
  gpm = GenotypePhenotypeMap(wildtype, genotypes, phenotypes,
                             stdeviations=stdeviations)

  # Show the dataFrame
  gpm.data

.. image:: _img/basic-example-df.png
    :width: 350px


Documentation
-------------

.. toctree::
   :maxdepth: 2

   _pages/simulate
   _pages/io
   _api/main.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
