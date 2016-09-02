
#Python API for analyzing and manipulating genotype-phenotype maps

This package defines a strict data-structure for genotype-phenotype mapping data.
It, then, provides a simple API for analyzing and manipulating such data. Some of
the things you can do:  

1. Convert any set of sequences into binary representations for modeling.
2. Paired with LatticeGPM, can easily construct lattice protein sequence spaces.
3. Seamlessly construct [NetworkX](https://networkx.github.io) object, enabling graph and network analysis.
4. Visualize networks with [NetworkViewer](https://github.com/harmslab/NetworkViewer) application.

## Installation

To install this package, clone from source and use pip.

```
git clone https://github.com/harmslab/seqspace
cd seqspace
pip install -e .
```

## Dependencies

The following modules are required for this to work. Also, the examples/tutorials are written in Jupyter notebooks and require IPython to be install.

* [NetworkX](https://networkx.github.io/)
* [Numpy](http://www.numpy.org/)
