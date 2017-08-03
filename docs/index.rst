.. gpmap documentation master file, created by
   sphinx-quickstart on Fri Jul  8 10:41:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

The **gpmap** package standardizes a data structure for genotype-phenotype (GP) maps.
Subset, manipulate, extend, etc. genotype-phenotype maps easily. Calculate statistics,
model evolutionary trajectories, predict phenotypes. Efficient memory usage and manipulation,
using Pandas Dataframe/Series.

This package includes modules for simulating computational genotype-phenotype maps
using methods described in the literature. See the Simulating_ page.

.. _Simulating: _pages/simulate.html

The GenotypePhenotypeMap object can be easily ported to network graphs (via NetworkX and GPGraph).

.. image:: _img/gpm.png
    :align: center

.. toctree::
   :maxdepth: 2

   _pages/gpm
   _pages/simulate
   _pages/io
   _api/main.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
