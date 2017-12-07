
# GPMap

*A Python API for managing genotype-phenotype map dat*

[![Join the chat at https://gitter.im/harmslab/gpmap](https://badges.gitter.im/harmslab/gpmap.svg)](https://gitter.im/harmslab/gpmap?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Documentation Status](https://readthedocs.org/projects/gpmap/badge/?version=latest)](http://gpmap.readthedocs.io/en/latest/?badge=latest)

Defines a flexible object for managing genotype-phenotype (GP) map data. At it's core,
`gpmap` stores the data in Pandas DataFrames and thus, interacts seamlessly with the
PyData egosystem.

<img src="docs/_img/gpm.png" align="middle">

## Basic example

Import the package's base object.
```python
from gpmap import GenotypePhenotypeMap
```

Load a dataset from disk.
```python
gpm = GenotypePhenotypeMap.read_json("data.json")
```

## Installation

To install this package, clone from source and use pip.
```
git clone https://github.com/harmslab/gpmap
cd gpmap
pip install -e .
```

## Dependencies

The following modules are required. Also, the examples/tutorials are written in Jupyter notebooks and require IPython to be install.

* [Numpy](http://www.numpy.org/)
* [Pandas](https://pandas.pydata.org/)
