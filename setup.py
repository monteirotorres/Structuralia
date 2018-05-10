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
  version = '0.1',
  description = 'A toolset to work with proteins',
  author = 'Pedro Torres',
  author_email = 'monteirotorres@gmail.com',
  url = 'https://github.com/monteirotorres/Structuralia', 
  download_url = 'https://github.com/monteirotorres/Structuralia/archive/0.1.tar.gz', 
  keywords = ['PDB', 'structure', 'protein'], 
  classifiers = [],
)
