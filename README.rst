Structuralia
############

This suite of python scripts was built to access the PDB database and both prepare inputs and read generated output files.

It relies mainly on Biopython and Biopandas to manipulate PDB structures in convenient ways.

Simply run AccessPDB to perform some of these functions intuitively.

The Toolbox is packed up with functions that can be imported and used in your own scripts.

Installation
************
The scripts are available as a `PyPi project`_. Just install them with:

.. _`PyPi project`: https://pypi.org/project/Structuralia/



:code:`pip install Structuralia`



Prerequisites
*************

The following packages and external programs are used by Structuralia scripts and must be installed and in either the binaries path or python path.

Python packages
===============

  - progressbar
  - pandas
  - biopython
  - biopandas
  - pathlib
  - parasail

External software (must be installed separately)
================================================

  - `TM-align`_
  - `GESAMT`_

.. _`TM-align`: https://zhanglab.ccmb.med.umich.edu/TM-align/
.. _`GESAMT`: http://www.ccp4.ac.uk/html/gesamt.html


Author
******

Pedro Torres, Ph.D

Department Of Biochemistry
University of Cambridge
80 Tennis Court Road
Cambridge CB2 1GA

License
*******

This project is licensed under GNU license, provided along with the package - see `LICENSE`_.

.. _LICENSE: https://github.com/monteirotorres/Structuralia/blob/master/LICENSE
