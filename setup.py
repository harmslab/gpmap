from setuptools import setup

setup(name='gpm',
      version='0.1',
      description='Internal mapping module for genotype-phenotype maps.',
      author='Zach Sailer',
      author_email='zachsailer@gmail.com',
      packages=['gpm'],
      install_requires=[
          'networkx',
          'numpy',
      ],
      zip_safe=False)