try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(name='seqspace',
      version='0.1',
      description='Internal mapping module for sequence space (specifically, genotype-phenotype maps).',
      author='Zach Sailer',
      author_email='zachsailer@gmail.com',
      packages=['seqspace'],
      install_requires=[
          'networkx',
          'numpy',
      ],
      zip_safe=False)