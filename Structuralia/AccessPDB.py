#!/opt/anaconda3/bin/python
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
# Imports
###############################################################################
import os
import textwrap as tw
import Structuralia.Structuralia.Toolbox as strtools


# Main Function
###############################################################################

def main():
    strtools.set_globals()
    print(tw.dedent("""\
          Structuralia  Copyright (C) 2018  Pedro H. M. Torres
          This program comes with ABSOLUTELY NO WARRANTY

          This program is intended to manipulate PDB files retreived from the
          PDB databank. Please run inside the directory where the input files
          are and where the output files and directories should be created.

          Valid options include:

          1) Download non-redundant pdb files (BlastClust).
          2) Extract a list of non-redundant PDBs from a given list.
          3) Download a list of pdb files.
          4) Clean the downloaded files and, select homo oligomers and sort
          them according to oligomeric state.
          5) Simply clean the PDB files in a directory.
          6) Select a single chain from PDB files in the directory and write
          them in a sub-directory

          """))
    pdb_dir = os.getcwd()
    option = input('Chose one of the options: ')
    if option == '1':
        cutoff = input('\nChose a cutoff. Available options are: 30, 40, 50, 70, 90, 95 and 100.\n')
        strtools.download_nr(pdb_dir, cutoff)
    elif option == '2':
        list_file = pdb_dir+'/'+input('Specify a file contining the list of pdbs: ')
        strtools.extract_nr(list_file)
    elif option == '3':
        list_file = pdb_dir+'/'+input('Specify a file contining the list of pdbs: ')
        strtools.download_pdblist(pdb_dir, list_file)
    elif option == '4':
        strtools.clean_and_sort(pdb_dir)
    elif option == '5':
        strtools.clean_pdb_files(pdb_dir)
    elif option == '6':
        strtools.single_chain(pdb_dir)
    else:
        print('\n Sorry, you did not select a valid option.\n\n Try again.')


# Execute
###############################################################################


main()
