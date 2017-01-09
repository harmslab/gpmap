try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(name='gpmap',
      version='0.1',
      description='Internal mapping module for sequence space (specifically, genotype-phenotype maps).',
      author='Zach Sailer',
      author_email='zachsailer@gmail.com',
      packages=['gpmap'],
      install_requires=[
          'networkx',
          'numpy',
      ],
      zip_safe=False)
