# License
###############################################################################
'''
Structuralia: A suite of python scripts to easily manipulate PDBs
Copyright (C) 2018  Pedro H. M. Torres

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact info:
Department Of Biochemistry
University of Cambridge
80 Tennis Court Road
Cambridge CB2 1GA
E-mail address: monteirotorres@gmail.com
'''

# Description
###############################################################################
'''
Set most widely used global variables across the Structuralia scripts.
'''

# Imports
###############################################################################
import os
import progressbar
import Bio.PDB as bpp

# Global Variables
###############################################################################

global workdir
global oligo_dict
global p
global io
global widgets
workdir = os.getcwd()+'/'
oligo_dict = {1: 'MONOMERIC', 2: 'DIMERIC', 3: 'TRIMERIC', 4: 'TETRAMERIC',
              5: 'PENTAMERIC', 6: 'HEXAMERIC'}
p = bpp.PDBParser(PERMISSIVE=0, QUIET=True)
io = bpp.PDBIO()
widgets = [' [', progressbar.SimpleProgress(), '] ',
           progressbar.Bar(),
           progressbar.Percentage(),
           ' (', progressbar.AdaptiveETA(), ') ']
