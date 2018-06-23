#!/usr/bin/env python3
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
This code contains several tools related to the determination of the oligomeric
state of complexes as well as tools to correlate the oligomeric state with
sequence identity.

Developed by: Pedro Torres
'''
# Imports
###############################################################################
import os
import re
import gzip
import Structuralia.Toolbox as strtools
from Structuralia.GlobalVars import *
import itertools as it
import Bio.pairwise2 as bpw2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from progressbar import progressbar as pg
import pandas as pd
import progressbar
import textwrap as tw
import subprocess
import numpy as np
import parasail
from functools import partial
from multiprocessing import Pool

# Dictionaries
###############################################################################

# Global Variables
###############################################################################
matrix = matlist.blosum62

# Classes
###############################################################################

# Functions
###############################################################################
'''
Functions to create alignment matrix in the form of a stack
(list of sublists) in which each line (sublist) contains the following
structure:

[id-1, id-2, percent-id, oligo_equivalence]

Note that oligo_equivalence is a binary variable which evaluates to
1 when the oligomeric state of both sequences is the same.
'''
def make_stack_simple(pair, id_seq_oligo_list):
    '''
    Using simple pairwise module from Biopython.
    Called by: OligoState.py:main()
    '''
    with open('stack.dat', 'a') as stack:
        a=pair[0]
        b=pair[1]
        line = []
        line.append(id_seq_oligo_list[a][0])
        line.append(id_seq_oligo_list[b][0])
        alignment = bpw2.align.globalxx(id_seq_oligo_list[b][1],
                                    id_seq_oligo_list[a][1],
                                    one_alignment_only=True)
        percent_id = alignment[0][2]*100/alignment[0][4]
        line.append(percent_id)
        if id_seq_oligo_list[a][2] == id_seq_oligo_list[b][2]:
            line.append(1)
        else:
            line.append(0)
        stack.writelines(str(line[0])+','+str(line[1])+','+str(line[2])+','+str(line[3])+'\n')


def make_stack_matrix(pair, id_seq_oligo_list):
    '''
    Using pairwise module from Biopython and employing blosum62
    substitution matrix.
    Called by: OligoState.py:main()
    '''
    with open('stack.dat', 'a') as stack:
        a=pair[0]
        b=pair[1]
        line = []
        line.append(id_seq_oligo_list[a][0])
        line.append(id_seq_oligo_list[b][0])
        alignment = bpw2.align.globalds(id_seq_oligo_list[b][1],
                                        id_seq_oligo_list[a][1],
                                        matrix,
                                        -1.0,
                                        -0.2,
                                        one_alignment_only=True)
        contents = format_alignment(*alignment[0])
        matches = len(re.findall('\|', contents))
        percent_id = matches*100/alignment[0][4]
        line.append(percent_id)
        if id_seq_oligo_list[a][2] == id_seq_oligo_list[b][2]:
            line.append(1)
        else:
            line.append(0)
        stack.writelines(str(line[0])+','+str(line[1])+','+str(line[2])+','+str(line[3])+'\n')


def make_stack_gesamt(pair, id_seq_oligo_list):
    '''
    Using structural alignment (gesamt) to calculate identities.
    Called by: OligoState.py:main()
    '''
    with open('stack.dat', 'a') as stack:
        a=pair[0]
        b=pair[1]
        apdb = str(id_seq_oligo_list[a][0])+'.pdb'
        bpdb = str(id_seq_oligo_list[b][0])+'.pdb'
        line = []
        line.append(id_seq_oligo_list[a][0])
        line.append(id_seq_oligo_list[b][0])
        try:
            percent_id = float(str(subprocess.check_output(['gesamt', apdb, bpdb])).split('Sequence Id:     : ')[1].split(' ')[0])
        except:
            percent_id = float(0)
        line.append(percent_id)
        if id_seq_oligo_list[a][2] == id_seq_oligo_list[b][2]:
            line.append(1)
        else:
            line.append(0)
        stack.writelines(str(line[0])+','+str(line[1])+','+str(line[2])+','+str(line[3])+'\n')


def make_stack_parasail(pair, id_seq_oligo_list):
    '''
    Using simple pairwise module from Biopython.
    Called by: OligoState.py:main()
    '''
    with open('stack.dat', 'a') as stack:
        a=pair[0]
        b=pair[1]
        line = []
        line.append(id_seq_oligo_list[a][0])
        line.append(id_seq_oligo_list[b][0])
        alignment = parasail.nw_stats_striped_16(id_seq_oligo_list[b][1],
                                                 id_seq_oligo_list[a][1],
                                                 10, 1,
                                                 parasail.blosum62)
        percent_id = (alignment.matches)/alignment.length*100
        line.append(percent_id)
        if id_seq_oligo_list[a][2] == id_seq_oligo_list[b][2]:
            line.append(1)
        else:
            line.append(0)
        stack.writelines(str(line[0])+','+str(line[1])+','+str(line[2])+','+str(line[3])+'\n')



def make_diagonal(id_seq_oligo_list, entry):
    '''
    Small function to add the diagonal values to the stack, wich are not
    computed by the above functions since itertools.combinations does not
    iterate through them.
    Called by: OligoState.py:main()
    '''
    line = []
    line.append(id_seq_oligo_list[entry][0])
    line.append(id_seq_oligo_list[entry][0])
    percent_id = 100
    line.append(percent_id)
    line.append(1)
    return line


def bin_stack(stackfile):
    binned_stacks = []
    for bin in range(1, 101):
        binned_stacks.append([line for line in open('stack.dat').readlines() if bin-1 < float(line.split(',')[2]) < bin])
    return binned_stacks

# Main Function
###############################################################################


def main():
    print(tw.dedent("""\
          Structuralia  Copyright (C) 2018  Pedro H. M. Torres
          This program comes with ABSOLUTELY NO WARRANTY

                               *** OligoState ***

          This tool iterates through the PDB files in current directory
          and correlates identity with oligomeric state.

          Chose a alignment scheme:

          1) Simple pairwise alignment, no gap penalties (faster)
          2) Employ the Blosum62 matrix (slower)
          3) Use Gesamt structural alignment to calculate identities (slowest)
          4) Use Parasail to calculate identities (fastest)
          """))
    option = input('Chose one of the options: ')
    if option != '1' and option !='2' and option != '3' and option != '4':
        print('\nNot a valid option.\n')
        exit()
    prefix = input('\nPlease enter an output files prefix: ')
    pdb_dir = os.getcwd()
    with open('stack.dat', 'w+'):
        pass
    id_seq_oligo_list = []
    for pdb in os.listdir():
        if pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz"):
            id_seq_oligo = []
            pdb_name = pdb.split('.')[0].split("/")[-1]
            if pdb.endswith(".ent.gz"):
                pdb_file = gzip.open(pdb_dir+'/'+pdb, 'rt')
            else:
                pdb_file = open(pdb_dir+'/'+pdb)
            try:
                structure = p.get_structure(pdb_name, pdb_file)
            except:
                print("Structure "+pdb_name+" could not be strictly parsed.")
                continue
            seq = strtools.extract_seqs(structure, 0)[1][0]
            if len(seq) > 69:
                id_seq_oligo.append(pdb_name)
                id_seq_oligo.append(seq)
                id_seq_oligo.append(pdb_name.split('-')[0])
                id_seq_oligo_list.append(id_seq_oligo)
    pairs = list(it.combinations(range(len(id_seq_oligo_list)), 2))
    po = Pool()
    if option == '1':
        f=partial(make_stack_simple, id_seq_oligo_list=id_seq_oligo_list)
        po.imap(f, pairs)
    elif option == '2':
        f=partial(make_stack_matrix, id_seq_oligo_list=id_seq_oligo_list)
        po.imap(f, pairs)
    elif option == '3':
        f=partial(make_stack_gesamt, id_seq_oligo_list=id_seq_oligo_list)
        po.imap(f, pairs)
    elif option == '4':
        f=partial(make_stack_parasail, id_seq_oligo_list=id_seq_oligo_list)
        po.imap(f, pairs)
    po.close()
    po.join()
    with open('stack.dat', 'a') as stack:
        for entry in range(len(id_seq_oligo_list)):
            line=make_diagonal(id_seq_oligo_list, entry)
            stack.writelines(str(line[0])+','+str(line[1])+','+str(line[2])+','+str(line[3])+'\n')
    stack = pd.read_csv('stack.dat', header=None)
    df = pd.DataFrame(stack, index=None, columns=None)
    df = pd.crosstab(index=df[0],columns=df[1], values=df[2], aggfunc='sum' , margins=False, dropna=True).fillna(0)
    df = np.around(df, decimals=2)
    df = df+df.T-np.diag(df.values.diagonal())
    df.to_csv(prefix+'matrix.csv')
    binned_stacks = bin_stack('stack.dat')
    nbin = 0
    nbins = []
    bin_counts = []
    same_state_list = []
    for bin in binned_stacks:
        nbin += 1
        nbins.append(nbin)
        bin_count = len(bin)
        bin_counts.append(bin_count)
        same_state = 0
        for line in bin:
            if line[3] == 1:
                same_state +=1
        if bin_count != 0:
            same_state_percent = same_state*100/bin_count
        else:
            same_state_percent = 0
        same_state_list.append(same_state_percent)
    final = pd.DataFrame({'Identity Range': nbins,
                          'Counts': bin_counts,
                          'Oligomerization Equivalence': same_state_list})
    final.to_csv(prefix+'.csv')


# Execute
###############################################################################
main()
