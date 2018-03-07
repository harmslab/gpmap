try:
    from setuptools import setup
except:
    from distutils.core import setup

# Packages and subpackages to install
packages = [
    'gpmap',
    'gpmap.simulate'
]

setup(name='gpmap',
      version='0.4.1',
      description='Data-structure for analyzing genotype-phenotype map data.',
      author='Zach Sailer',
      author_email='zachsailer@gmail.com',
      packages=packages,
      url="https://github.com/harmslab/gpmap",
      install_requires=[
          'numpy',
          'scipy',
          'pandas'
      ],
      classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
        ],
      zip_safe=False)
