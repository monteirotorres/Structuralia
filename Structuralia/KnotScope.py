#!/usr/bin/env python3
# Imports
###############################################################################
import os
import sys
import time
import argparse
import collections
import pandas as pd
import textwrap as tw
import itertools as it
import Structuralia.Toolbox as strtools
from Structuralia.GlobalVars import *
import Bio.PDB.Polypeptide as bpp_poly

# LICENSE
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
 A simple biopython-based code to detect and output unrealistic aplha
 carbon distances in pdb files.
'''
# Dictionaries
###############################################################################

# Global Variables
###############################################################################

# Classes
###############################################################################

# Functions
###############################################################################
def printv(text, verbosity):
    '''
    Function to control print verbosity.
    Called by: get_knots()
    '''
    if verbosity == 1:
        print(text)


def make_diagonal(core, entry):
    '''
    Outputs a list of 3 elements to act as the diagonal of an identity matrix.
    the first and two elements are equal and will serve as the row and column
    indexes, while the third element is a dash.
    Called by: get_knots()
    '''
    line = []
    line.append(core[entry].get_parent().id[1])
    line.append(core[entry].get_parent().id[1])
    line.append('-')
    return line


def get_knots(pdb, cutoff, cluster_cutoff, genpdb, verbosity):
    '''
    Main routine, uses biopython and pandas to detect knots and cluster them
    through the implementation of the average linkage algorithm
    Called by: main()
    '''
    if (pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz") or pdb.endswith(".pdb1") or pdb.endswith(".pdb1.gz") or pdb.endswith(".pdb.gz")) and not pdb.startswith('CONTACTS-'):
        pdb_name, structure, nchains = strtools.parse_pdb_structure(pdb)
        print(str('\n'+clrs['p']+pdb+clrs['n']))
        with open('KnotScope.log', 'a') as log:
            log.write(str('\n[STRUCTURE],'+pdb+'\n'))
        mainchain = [atom for atom in bpp.Selection.unfold_entities(structure[0], 'A') if bpp_poly.is_aa(atom.get_parent(), standard=True) and (atom.id == 'CA')]# or atom.id == 'N' or atom.id == 'O')]
        contacts = []
        core = []
        for atom in mainchain:
            distances = []
            ns = bpp.NeighborSearch(mainchain)
            center = atom.get_coord()
            neighbors = [neighbor for neighbor in ns.search(center, cutoff) if (neighbor.get_parent().id[1] - atom.get_parent().id[1]) > abs(3)]
            if neighbors:
                for neighbor in neighbors:
                    d = neighbor - atom
                    distances.append(d)
                    if d <= cutoff:
                        printv(clrs['y']+'Unlikely proximity'+clrs['n']+' between residues '+clrs['y']+str(atom.get_parent().id[1])+clrs['n']+' and '+clrs['y']+str(neighbor.get_parent().id[1])+clrs['n']+'!', verbosity)
                        printv(str(d), verbosity)
                        with open('KnotScope.log', 'a') as log:
                            log.write('[CLASH],'+str(atom.get_parent().id[1])+','+str(neighbor.get_parent().id[1])+','+str(d)+'\n')
                        contacts.append(neighbor.get_parent())
                        contacts.append(atom.get_parent())
                        if atom not in core:
                            core.append(atom)
                        if neighbor not in core:
                            core.append(neighbor)
        # Save contacts to pdb file if they exist
        if contacts and genpdb:
            io.set_structure(structure)
            io.save('CONTACTS-'+pdb, strtools.SelectResidues(contacts))
        # Start cluster analysis to separate knots
        pairwisedist = []
        # Measure pairwise distances of every CA involved in knots and record in vertical list
        if len(core) > 1:
            for a, b in it.combinations(core, 2):
                d = a - b
                pairwisedist.append([a.get_parent().id[1],b.get_parent().id[1],d])
            # Add values for diagonal
            for entry in range(len(core)):
                line=make_diagonal(core, entry)
                pairwisedist.append([line[0],line[1],line[2]])
            # Create pandas dataframe, make it a square and symmetric matrix
            df = pd.DataFrame(pairwisedist, index=None, columns=None)
            df = pd.crosstab(index=df[0],columns=df[1], values=df[2], aggfunc='sum' , dropna=True).fillna(0)
            df = df+df.T
            # Start average linkage algorithm
            reslist = list(df.columns)
            clusters = []
            row_index = -1
            col_index = -1
            array = []
            for n in range(df.shape[0]):
                array.append(n)
            clusters.append(array.copy())
            for k in range(1, df.shape[0]):
                min_val = sys.maxsize
                for i in range(0, df.shape[0]):
                    for j in range(0, df.shape[1]):
                        #print(str(df.iloc[i,j]))
                        if type(df.iloc[i, j]) != str:
                            if(df.iloc[i, j] <= min_val):
                                min_val = df.iloc[i,j]
                                row_index = i
                                col_index = j

                for i in range(0,df.shape[0]):
                    if(i != col_index and i!=row_index):
                        temp = (df.iloc[col_index,i]+df.iloc[row_index,i])/2
                        df.iloc[col_index,i] = temp
                        df.iloc[i,col_index] = temp
                for i in range (0,df.shape[0]):
                    df.iloc[row_index,i] = sys.maxsize
                    df.iloc[i,row_index] = sys.maxsize
                minimum = min(row_index,col_index)
                maximum = max(row_index,col_index)
                for n in range(len(array)):
                    if(array[n] == maximum):
                        array[n] = minimum
                clusters.append(array.copy())
                # Stop iterations when minimum pairwise distance in the matrix is greater than 22
                if min_val > cluster_cutoff:
                    break
            # Get the clusters from last iteration and 'count' elements
            clustered_res = clusters[-1]
            counter = collections.Counter(clustered_res)
            # Combine residue and cluster information and print them user-friendly
            clusterdict = dict(zip(reslist, clustered_res))
            print(clrs['y']+'\nLikely '+str(len(set(clusters[-1])))+' knot(s) found in structure under chosen criteria...'+clrs['n'])
            n = 0
            k_lengths = []
            for cl in set(clusterdict.values()):
                n += 1
                cluster_residues = []
                print('\nKnot '+clrs['y']+str(n)+clrs['n']+' (Cluster id: '+str(cl)+') involves '+clrs['y']+str(list(counter.values())[n-1])+clrs['n']+' residues:')
                for res in clusterdict:
                    if clusterdict[res] == cl:
                        cluster_residues.append(res)
                k_lengths.append(len(cluster_residues))
                print(clrs['y']+', '.join([str(a) for a in cluster_residues])+clrs['n'])
                with open('KnotScope.log', 'a') as log:
                    log.write('[K'+str(n)+'-RES],'+','.join([str(a) for a in cluster_residues])+'\n')
                    log.write('[K'+str(n)+'-LEN],'+str(len(cluster_residues))+'\n')
            with open('KnotScope.log', 'a') as log:
                log.write('[SUM],str,'+pdb+',ca_clash,'+str(len(core))+',nknots,'+str(len(set(clusters[-1])))+',longest,'+str(max(k_lengths))+'\n')
            return clusterdict
        else:
            print(clrs['g']+'No CA distances under '+str(cutoff)+' angstrons found'+clrs['n']+'!\n')
            with open('KnotScope.log', 'a') as log:
                log.write('[SUM],str,'+pdb+',ca_clash,0,nknots,0,longest,0\n')

        del pdb_name, structure, nchains, contacts
    elif pdb.startswith('CONTACTS-'):
        pass
    else:
        print(clrs['y']+pdb+clrs['n']+' not a pdb-related structure format.'+clrs['r']+' SKIPPING!'+clrs['n'])


# Main Function
###############################################################################
def main():
    print(tw.dedent("""
         ###################################################################
         ########################### KnotScope #############################
         ###################################################################

                 Copyright (C) 2018  Torres, P.H.M.; Blundell, T.L.
                          [The University of Cambridge]

                   This program comes with ABSOLUTELY NO WARRANTY

         A simple biopython-based code to detect and output unrealistic aplha
         carbon distances in pdb files.

         ###################################################################

          """))

    parser = argparse.ArgumentParser()

    parser.add_argument(dest='files',
                        metavar='files',
                        help='Either pdb file or directory containing pdbs\n')

    parser.add_argument('-v', '--verbose',
                        dest='verbosity',
                        action='count',
                        default=0,
                        help='Controls verbosity\n')

    parser.add_argument('-p', '--genpdb',
                        dest='genpdb',
                        action='count',
                        default=0,
                        help='Defines whether pdbs should be generated.\n')

    parser.add_argument('-c', '--cutoff',
                        dest='cutoff',
                        type=float, default=3.2,
                        metavar='float',
                        help='Cutoff to be used in neighborsearching (default = 3.2)\n')

    parser.add_argument('-l', '--cluster',
                        dest='cluster_cutoff',
                        type=float, default=20.0,
                        metavar='float',
                        help='Cutoff to be used in knot clustering (default = 20)\n')

    args = parser.parse_args()

    files = args.files
    verbosity = args.verbosity
    cutoff = args.cutoff
    cluster_cutoff = args.cluster_cutoff
    genpdb = args.genpdb
    if genpdb == 1:
        genpdb = True
    elif genpdb == 0:
        genpdb = False

    now = str(time.strftime("%d-%m-%Y@%H.%M.%S"))

    with open('KnotScope.log', 'w+') as f:
        f.write('[KNOTSCOPE LOG],'+now+'\n')
        f.write('[CRITERIA],cutoff,'+str(cutoff)+',cluster_cutoff,'+str(cluster_cutoff)+',genpdb,'+str(genpdb)+'\n')
    if os.path.isdir(args.files):
        for pdb_file in os.listdir(args.files):
            get_knots(pdb_file, cutoff, cluster_cutoff, genpdb, verbosity)
    elif os.path.isfile(args.files):
        get_knots(args.files, cutoff, cluster_cutoff, genpdb, verbosity)
    else:
        raise FormatError('Please provide either a pdb file or a valid path to files')


# Execute
###############################################################################
if __name__ == "__main__":
    main()
