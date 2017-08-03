.. gpmap documentation master file, created by
   sphinx-quickstart on Fri Jul  8 10:41:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

The `gpmap` package standardizes a data structure for genotype-phenotype (GP) maps.
Subset, manipulate, extend, etc. genotype-phenotype maps easily. Calculate statistics,
model evolutionary trajectories, predict phenotypes. Efficient memory usage and manipulation,
using Pandas Dataframe/Series.

This package includes methods to simulate various computatioanl genotype-phenotype maps
present in the literature.

The GenotypePhenotypeMap object can is easily ported to network graphs (via NetworkX and GPGraph).

.. image:: _img/gpm.png
    :align: center

.. toctree::
   :maxdepth: 2

   _pages/gpm
   _pages/io
   _pages/tutorials
   _pages/simulate
   _api/main.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
