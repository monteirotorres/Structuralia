from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name = 'Structuralia',
  packages = find_packages(),
  version = '0.9',
  description = 'A toolset to work with proteins',
  long_description = long_description,
  author = 'Pedro Torres',
  author_email = 'monteirotorres@gmail.com',
  url = 'https://github.com/monteirotorres/Structuralia',
  download_url = 'https://github.com/monteirotorres/Structuralia/archive/0.9.tar.gz',
  keywords = ['PDB', 'structure', 'protein'],
  classifiers = [],
  py_modules=[
      'AccessPDB',
      'OligoSum',
      'OligoState'
      'Toolbox'],
  entry_points={
      'console_scripts': [
          'AccessPDB = Structuralia.AccessPDB:main',
          'OligoSum = Structuralia.OligoSum:main',
          'OligoState = Structuralia.OligoState:main']},
  install_requires=[
      'markdown',
      'biopandas',
      'pandas',
      'progressbar2',
      'pathlib2',
      'biopython'])
