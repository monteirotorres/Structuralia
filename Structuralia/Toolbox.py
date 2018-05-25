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
This file contains several functions and is supposed to be imported by other
python programs. Several other Structuralia programs will rely upon it.
The recommended import command is:

import Structuralia.Toolbox as strtools

Developed by: Pedro Torres
'''

# Imports
###############################################################################
import os
import re
import gzip
import urllib
import shutil
import subprocess
import progressbar
import pandas as pd
import textwrap as tw
import itertools as it
import Bio.PDB as bpp
import Bio.pairwise2 as bpw2
import Bio.PDB.Polypeptide as bpp_poly
from biopandas.pdb import PandasPdb
from progressbar import progressbar as pg
from Structuralia.TimeWrap import timed
from Structuralia.GlobalVars import *


# Classes
###############################################################################


class SelectAA(bpp.Select):
    '''
    Biopython select class to select only aminoacids
    Called by: clean_pdb()
    '''
    def accept_residue(self, residue):
        if bpp_poly.is_aa(residue.get_resname(), standard=True):
            return 1
        else:
            return 0


class SelectChain(bpp.Select):
    '''
    Biopython select class to select a single chain id.
    Called by: single_chain()
    '''
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        if chain.id == self.chain_id:
            return 1
        else:
            return 0


# Functions
###############################################################################


# @timed
def get_nr_list(cutoff):
    '''
    Uses Biopython to retrieve a list of nonredundant PDBs from RCSB
    Called by: download_nr30()
               extract_nr()
    '''
    url = 'ftp://resources.rcsb.org/sequence/clusters/bc-'+cutoff+'.out'
    pdblist = []
    for entry in urllib.request.urlopen(url):
        pdblist.append(str(entry.decode('utf-8')).lstrip().split('_')[0])
    nr_list = list(set(pdblist))
    return nr_list


# @timed
def read_list_file(list_file):
    '''
    Simply reads a text file with a single pdb code per line and
    makes a list with them (changing lower to uppercase, if necessary).
    Called by: download_pdblist()
               extract_nr()
    '''
    pdb_list = []
    with open(list_file, 'r') as f:
        pdb_list = f.read().splitlines()
    pdb_list_upper = [i.upper() for i in pdb_list]
    return pdb_list_upper


# @timed
def download_nr(pdb_dir, cutoff):
    '''
    Downloads representative PDB files into
    specified directory. (Option 1)
    Called by: AccessPDB.py:main()
    '''
    pdb_codes = get_nr_list(cutoff)
    pdbl = bpp.PDBList()
    pdbl.download_pdb_files(pdb_codes, pdir=pdb_dir, file_format='pdb')


# @timed
def extract_nr(list_file):
    '''
    Extract a list of non-redundant PDBs from a given list. (Option 2)
    Called by: AccessPDB.py:main()
    '''
    nr_list = get_nr_list()
    pdb_list = read_list_file(list_file)
    new_nr_list = list(set(nr_list).intersection(pdb_list))
    new_list_file = list_file+".nr"
    with open(new_list_file, 'a') as f:
        for entry in new_nr_list:
            f.write(entry+'\n')


# @timed
def download_pdblist(pdb_dir, list_file):
    '''
    Uses Biopython to download a list of PDBs into specified directory.
    (Option 3)
    Called by: AccessPDB.py:main()
    '''
    pdb_codes = read_list_file(list_file)
    pdbl = bpp.PDBList()
    pdbl.download_pdb_files(pdb_codes, pdir=pdb_dir, file_format='pdb')

def parse_pdb(pdb_dir, pdb):
    if pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz"):
        pdb_name = pdb.split('.')[0].split("/")[-1]
        if pdb.endswith(".ent.gz"):
            pdb_file = gzip.open(pdb_dir+'/'+pdb, 'rt')
            contents = gzip.open(pdb_dir+'/'+pdb, 'rt').read()
        else:
            pdb_file = open(pdb_dir+'/'+pdb)
            contents = open(pdb_dir+'/'+pdb, 'rt').read()
        try:
            structure = p.get_structure(pdb_name, pdb_file)
        except:
            print("Structure "+pdb_name+" could not be strictly parsed.")
    return pdb_name, pdb_file, structure, contents

# @timed
def extract_seqs(structure, defmodel):
    '''
    Uses Biopython to count the numer of chains and to extract the
    each chain's sequence as a list of sequences.
    Called by: clean_and_sort()
    '''
    nchains = 0
    for model in structure:
        if model.id == defmodel:
            seqs = []
            chain_ids = []
            for chain in model:
                nchains += 1
                seqlist = []
                for residue in chain:
                    if bpp_poly.is_aa(residue.get_resname(),
                                      standard=True):
                        seqlist.append(
                            bpp_poly.three_to_one(residue.get_resname()))
                    else:
                        seqlist.append('X')
                seq = str("".join(seqlist))
                seqs.append(seq)
                chain_ids.append(chain.id)
    return nchains, seqs, chain_ids


# @timed
def author_agrees(oligo_dict, contents, nchains):
    '''
    Searches in the original PDB file for the oligomeric
    status determined by the author.
    Called by: clean_and_sort()
    '''
    pattern = r"AUTHOR DETERMINED BIOLOGICAL UNIT: "+oligo_dict[nchains]
    if re.search(pattern, contents):
        return True
    else:
        return False


# @timed
def get_ids(seqs, nchains):
    '''
    Receives a list of sequences and calculates the identity matrix
    among them, as a list, using Biopython and itertools.
    Called by: clean_and_sort()
    '''
    ids = []
    for a, b in it.combinations(range(nchains), 2):
        score = bpw2.align.globalxx(seqs[a], seqs[b])[0][2]
        percent_id = score*100/len(seqs[a])
        ids.append(percent_id)
    return ids


# @timed
def clean_pdb(structure, pdb_name, clean_dir):
    '''
    Function to select and write pdb with only aminoacids
    Invokes SelectAA class constructed with Bio.PDB.select
    Called by: clean_pdb_files()
               clean_and_sort()
    '''
    reslist = []
    clean_name = clean_dir+pdb_name+'.clean.pdb'
    for res in structure.get_residues():
        if bpp_poly.is_aa(res.get_resname(), standard=True):
            reslist.append(res.resname)
    if len(reslist) > 30:
        io.set_structure(structure)
        io.save(clean_name, SelectAA())
        return True
    else:
        return False


# @timed
def clean_pdb_files(pdb_dir):
    '''
    Iterates over files in given directory and uses clean_pdb
    function to write PDB files containing only amino_acids. (Option 5)
    Called by: AccessPDB.py:main()
    '''
    clean_dir = "clean"
    os.mkdir(clean_dir)
    print('\n\nCleaning PDB files...\n')
    for pdb in pg(os.listdir(pdb_dir), widgets=widgets):
        if pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz"):
            pdb_name = pdb.split('.')[0].split("/")[-1]
            if pdb.endswith(".ent.gz"):
                pdb_file = gzip.open(pdb_dir+'/'+pdb, 'rt')
            else:
                pdb_file = pdb_dir+'/'+pdb
            try:
                structure = p.get_structure(pdb_name, pdb_file)
            except:
                print("Structure "+pdb_name+" could not be strictly parsed.")
                continue
            clean_pdb(structure, pdb_name, clean_dir)


# @timed
def clean_and_sort(pdb_dir):
    '''
    Make clean directory and homo multimer subdirectories (Option 4)
    Called by: AccessPDB.py:main()
    '''
    clean_dir = "clean/"
    try:
        os.mkdir(clean_dir)
        for i in range(1, 7):
            os.mkdir(clean_dir+'/'+str(i)+'mers')
    except:
        pass
    '''
    Loop through pdb files to detect homo get_oligomeric_status
    '''
    for pdb in pg(os.listdir(pdb_dir), widgets=widgets):
        if pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz"):
            pdb_name = pdb.split('.')[0].split("/")[-1]
            if pdb.endswith(".ent.gz"):
                pdb_file = gzip.open(pdb_dir+'/'+pdb, 'rt')
                contents = gzip.open(pdb_dir+'/'+pdb, 'rt').read()
            else:
                pdb_file = open(pdb_dir+'/'+pdb)
                contents = open(pdb_dir+'/'+pdb, 'rt').read()
            try:
                structure = p.get_structure(pdb_name, pdb_file)
            except:
                print("Structure "+pdb_name+" could not be strictly parsed.")
                continue
            nchains, seqs = extract_seqs(structure, 0)
            print("\n\nAssessing "+pdb_name+". This PDB has got "+str(nchains)+" chain(s).")
            if 2 <= nchains <= 6:
                if author_agrees(oligo_dict, contents, nchains):
                    print("Author agrees that "+pdb_name+" is "+oligo_dict[nchains]+" and IDs will be checked.")
                    ids = get_ids(seqs, nchains)
                    if all(id > 90 for id in ids):
                        print("All identities over 90%. Likely homo-oligomer. Cleaning and sorting.\n\n")
                        if clean_pdb(structure, pdb_name, clean_dir):
                            os.rename(clean_dir+pdb_name+'.clean.pdb',
                                      clean_dir+str(nchains)+'mers/'+pdb_name+'.clean.pdb')
                        else:
                            print("Oops! Polypeptide chain too short or inexistent. Skipping.\n\n")
                    else:
                        print("Identity under 90%. Likely not a homo-oligomer. Skipping.\n\n")
                else:
                    print("Author disagrees. Although PDB has "+str(nchains)+" chains, likely not "+oligo_dict[nchains]+".\n\n")
            elif nchains == 1:
                if author_agrees(oligo_dict, contents, nchains):
                    print("Author agrees that "+pdb_name+" is "+oligo_dict[nchains]+". Cleaning and sorting.\n\n")
                    if clean_pdb(structure, pdb_name, clean_dir):
                        os.rename(clean_dir+pdb_name+'.clean.pdb',
                                  clean_dir+str(nchains)+'mers/'+pdb_name+'.clean.pdb')
                    else:
                        print("Oops! Polypeptide chain too short or inexistent. Skipping.\n\n")
                else:
                    print("Author disagrees. Although PDB has "+str(nchains)+" chains, likely not "+oligo_dict[nchains]+".\n\n")
            elif nchains > 6:
                print("Too many chains. Skipping\n\n")


# @timed
def single_chain(pdb_dir):
    '''
    Iterates through a directory and uses Biopython to
    select and write the first chain from each pdb.
    Called by: AccessPDB.py:main() (Option 6)
    '''
    single_chain_dir = pdb_dir+"/SingleChains/"
    os.makedirs(single_chain_dir)
    print('\n\nExtracting first chain of each PDB file...\n')
    for pdb in pg(os.listdir(pdb_dir), widgets=widgets):
        if pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz"):
            pdb_name = pdb.split('.')[0].split("/")[-1]
            if pdb.endswith(".ent.gz"):
                pdb_file = gzip.open(pdb_dir+'/'+pdb, 'rt')
            else:
                pdb_file = pdb_dir+'/'+pdb
            try:
                structure = p.get_structure(pdb_name, pdb_file)
            except:
                print("Structure "+pdb_name+" could not be strictly parsed.")
                continue
            chains = structure.get_chains()
            for chain in chains:
                chain_id = chain.id
                break
            single_chain_name = single_chain_dir+pdb_name+chain_id+".pdb"
            io.set_structure(structure)
            io.save(single_chain_name, SelectChain(chain_id))


# @timed
def count_chains(structure):
    '''
    Uses Biopython to obtain the number of chains.
    Called by: OligoSum.py:main()
    '''
    nchains = 0
    for model in structure:
        if model.id == 0:
            for chain in model:
                nres = 0
                nchains += 1
                for residue in chain:
                    nres += 1
    return nchains, nres


# @timed
def merge_chains(pdb_file):
    '''
    Uses Biopandas to merge the chains of a pdb file into a single chain "A",
    renumber the residues, removes TER entries and saves a pdb file.
    Also returns the name of the created file.
    Called by: OligoSum.py:main()
    '''
    merged_file = pdb_file[:-4]+'.merged.pdb'
    prot = PandasPdb().read_pdb(pdb_file)
    chains = prot.df['ATOM']['chain_id'].unique()
    nres_list = []
    add_res = 0
    for chain in chains:
        slice = prot.df['ATOM'][prot.df['ATOM']['chain_id'] == chain]
        protdf = pd.DataFrame(prot.df['ATOM'])
        protdf.loc[protdf['chain_id'] == chain, ['residue_number']] = pd.to_numeric(protdf.loc[:,'residue_number'] + add_res, downcast='integer')
        nres = len(slice['residue_number'].unique())
        nres_list.append(nres)
        add_res = sum(nres_list)
    protdf['chain_id'] = 'A'
    prot._df['ATOM'] = protdf
    rm_ter = prot.df['OTHERS'][prot.df['OTHERS']['record_name'] != 'TER']
    prot._df['OTHERS'] = rm_ter
    if os.path.isfile(merged_file):
        print('\nReplacing old pdb file '+merged_file)
        os.remove(merged_file)
    prot.to_pdb(path=merged_file, records=None, gz=False, append_newline=True)
    return merged_file


# @timed
def run_tmalign(modeled, original):
    '''
    Externally runs the TMalign executable and returns the number of
    aligned residues, the RMSD and TM-Score values normalized by the
    length of the modeled structure.
    Called by: OligoSum.py:main()
    '''
    tm_result = str(subprocess.check_output(['TMalign', modeled, original]))
    tm_split = re.split(',|=|\n|\\\\n|\(', tm_result)
    tm_split
    aligned = int(tm_split[23])
    rmsd = float(tm_split[25])
    tmscore = float(tm_split[33])
    return aligned, rmsd, tmscore


def clean_and_sort_PDB(pdb_base):
    '''
    Runs the clean and sort function for the whole pdb database, considering
    the divided scheme adopted by RCSB, in which the subdirectories are
    the two middle characters in the PDB code.
    Called by: AccessPDB.py:main() (option 7)
    '''
    os.chdir(os.getcwd())
    for dir in os.listdir(pdb_base):
        subfolder = pdb_base+'/'+dir
        clean_and_sort(subfolder)

def min_chain_length(pdb_dir, length):
    '''
    Iterates over PDB files in directory, checks the chain lenght for each
    chain and copies the ones in which at least one chain has more than the
    desired length.
    Called by: AccessPDB.py:main() (option 8)
    '''
    filtered_dir = 'Over'+length
    os.mkdir(filtered_dir)
    for pdb in pg(os.listdir(pdb_dir), widgets=widgets):
        try:
            pdb_name, pdb_file, structure, contents = parse_pdb(pdb_dir, pdb)
        except:
            continue
        nchains, seqs = extract_seqs(structure, 0)
        if all(len(seq) < int(length) for seq in seqs):
            continue
        else:
            shutil.copyfile(pdb_dir+'/'+pdb_name+'.pdb',
                            pdb_dir+'/'+filtered_dir+'/'+pdb_name+'.pdb')

def pdb_to_fasta(pdb_dir):
    '''
    Simply iterates over the PDB files in the current directory and creates
    a FASTA file containing entries for each chain of all PDBs.
    Called by: AccessPDB.py:main() (option 9)
    '''
    fasta_file = pdb_dir.split("/")[-1]+'.fasta'
    for pdb in pg(os.listdir(pdb_dir), widgets=widgets):
        try:
            pdb_name, pdb_file, structure, contents = parse_pdb(pdb_dir, pdb)
        except:
            continue
        nchains, seqs, chain_ids = extract_seqs(structure, 0)
        with open(fasta_file, 'a') as f:
            for seq, chain_id in zip(seqs, chain_ids):
                wrapped_seq = "\n".join(tw.wrap(seq))
                fasta_entry = '>'+pdb_name+':'+str(chain_id)+'\n'+wrapped_seq+'\n\n'
                f.write(fasta_entry)
