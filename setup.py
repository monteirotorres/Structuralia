from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name = 'Structuralia',
  packages = ['Structuralia'],
  version = '0.5',
  description = 'A toolset to work with proteins',
  long_description = long_description,
  author = 'Pedro Torres',
  author_email = 'monteirotorres@gmail.com',
  url = 'https://github.com/monteirotorres/Structuralia',
  download_url = 'https://github.com/monteirotorres/Structuralia/archive/0.5.tar.gz',
  keywords = ['PDB', 'structure', 'protein'],
  classifiers = [],
  entry_points={
      'console_scripts': [
          'AccessPDB = AccessPDB:main',
          'OligoSum = OligoSum:main',
          'OligoState = OligoState:main']},
  install_requires=[
      'markdown',
      'biopandas',
      'pandas',
      'progressbar2',
      'pathlib2',
      'biopython'])
