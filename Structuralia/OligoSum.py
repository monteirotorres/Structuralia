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
This script was constructed to benchmark the oligomeric modelling pipeline,
previously developed by Sony Malhotra at the Biochemistry Department of the
Uiversity of Cambridge.

It must be run inside the top output directory and relies on a specific naming
scheme for the original complexes, to which the results shall be compared
(since it is a benchmark). The original complexes must be located one directory
up and have the ".clean." tag before the pdb file extension.

Developed by: Pedro Torres
'''

# Imports
###############################################################################
import re
import os
import pandas as pd
from progressbar import progressbar as pg
from pathlib import Path
import Structuralia.Toolbox as strtools
from Structuralia.GlobalVars import *


# Dictionaries
###############################################################################

# Global Variables
###############################################################################

benchmark_dir = '/home/torres/work/oligo/benchmarks/'

available_models_file = benchmark_dir+'available_models.lst'

clean_dir = str(Path(workdir).parents[0])+'/'

with open(available_models_file, 'r') as f:
    available_models = f.read().splitlines()

solutions_list = []

# Classes
###############################################################################

# Functions
###############################################################################

# Main Function
###############################################################################


def main():
    assert not os.path.isfile('OligoSum.csv'), '\033[1;31;40m \n\n File OligoSum.csv exists. Get rid of it.\n'
    results = open('OligoSum.csv', 'a')
    results.write('PDB ID,Chain length,Was Available,No of templates,Template,Is Same,ID,Gesamt Rank,Model Chains, Orig Chains, RMSD, Aligned (%), TM-Score\n')
    for job in pg(os.listdir(workdir), widgets=widgets):
        if job.endswith('homo.oligo'):
            solution_list = []
            tm_list = []
            ntemplates = 0
            template = 'NA'
            is_same = 'NA'
            id = 'NA'
            gesamt_rank = 'NA'
            RMSD = 'NA'
            alignedp = 'NA'
            nchains_model = 0
            job_id = job.split('_')[0]
            pdb_id = ''.join(list(job)[3:7])
            original_pdb_file = clean_dir+'pdb'+pdb_id+'.clean.pdb'
            original_structure = p.get_structure('original', original_pdb_file)
            nchains_orig, nres = strtools.count_chains(original_structure)
            if pdb_id in available_models:
                was_available = 'YES'
            else:
                was_available = 'NO'
            model_list = []
            for model in os.listdir(workdir+job):
                if model.startswith('oligo_model'):
                    model_list.append(model)
            if not model_list:
                tmscore = 'No templates found'
                tm_list.append(tmscore)
                results.write(pdb_id+
                              ','+str(nres)+
                              ','+was_available+
                              ','+str(ntemplates)+
                              ','+template+
                              ','+is_same+
                              ','+str(id)+
                              ','+str(gesamt_rank)+
                              ','+str(nchains_model)+
                              ','+str(nchains_orig)+
                              ','+str(RMSD)+
                              ','+str(alignedp)+
                              ','+str(tmscore)+'\n')
            else:
                for model in model_list:
                    template = model.split('_')[2]
                    if template == pdb_id:
                        is_same = 'YES'
                    else:
                        is_same = 'NO'
                    ntemplates = len(model_list)
                    gesamt_results_file = workdir+job+'/'+job_id+'.pdb_first_ges.res'
                    with open(gesamt_results_file, 'r') as f:
                        for line in f.readlines():
                            if re.search(template, line):
                                gesamt_rank = line.split()[0]
                                id = line.split()[4]
                                break
                    model_pdb_file = workdir+job+'/'+model+'/'+job_id+'.B99990001.pdb'
                    try:
                        model_structure = p.get_structure('modeled', model_pdb_file)
                        nchains_model = strtools.count_chains(model_structure)[0]
                    except FileNotFoundError:
                        nchains_model = 0
                    if nchains_model == nchains_orig:
                        merged_model_file = strtools.merge_chains(model_pdb_file)
                        merged_original_file = strtools.merge_chains(original_pdb_file)
                        aligned, RMSD, tmscore = strtools.run_tmalign(merged_model_file, merged_original_file)
                        alignedp = (aligned*100)/(nchains_model*nres)
                        tm_list.append(tmscore)
                    elif nchains_model == 0:
                        tmscore = 'Templates found but models not built'
                        RMSD = 'NA'
                        alignedp = 'NA'
                        tm_list.append(tmscore)
                    else:
                        tmscore = 'Wrong number of chains'
                        RMSD = 'NA'
                        alignedp = 'NA'
                        tm_list.append(tmscore)
                    results.write(pdb_id+
                                  ','+str(nres)+
                                  ','+was_available+
                                  ','+str(ntemplates)+
                                  ','+template+
                                  ','+is_same+
                                  ','+str(id)+
                                  ','+str(gesamt_rank)+
                                  ','+str(nchains_model)+
                                  ','+str(nchains_orig)+
                                  ','+str(RMSD)+
                                  ','+str(alignedp)+
                                  ','+str(tmscore)+'\n')

            for i in tm_list:
                    if type(i) == str or i != max([x for x in tm_list if type(x) != str]) or i < 0.58:
                        solution_list.append('NO')
                    else:
                        solution_list.append('YES')
            solutions_list.append(solution_list)
    results.close()
    merged_solutions_list = sum(solutions_list, [])
    results = pd.read_csv('OligoSum.csv', na_filter=False)
    results.insert(13, 'Solution', merged_solutions_list)
    results.to_csv('OligoSum.csv')


# Execute
###############################################################################

main()
