#!/usr/bin/env python
# File created on 30 Dec 2009.
from __future__ import division
from subprocess import call, check_call, CalledProcessError
from os import makedirs
from os.path import split, splitext

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

"""
This file contains the QIIME workflow functions which string together 
independent scripts. For usage examples see the related files in the 
scripts directory:
 - 
"""

## Begin functions used generally by the workflow functions
    
def print_commands(commands,status_update_callback):
    """Print list of commands to run """
    for c in commands:
        for e in c:
            status_update_callback('#%s' % e[0])
            print '%s' % e[1]
            
def call_commands_serially(commands,status_update_callback):
    """Run list of commands, one after another """
    for c in commands:
        for e in c:
            status_update_callback('%s\n%s' % e)
            try:
                check_call(e[1].split())
            except CalledProcessError, err:
                msg = "\n\n*** ERROR RAISED DURING STEP: %s\n" % e[0] +\
                 "Command run was:\n %s\n" % e[1] +\
                 "Command returned exit status: %d" % err.returncode
                print msg
                exit(-1)

def print_to_stdout(s):
    print s
    
def no_status_updates(s):
    pass

def get_params_str(params):
    result = []
    for param_id, param_value in params.items():
        result.append('--%s' % (param_id))
        if param_value != None:
            result.append(param_value)
    return ' '.join(result)

## End functions used generally by the workflow functions

## Begin task-specific workflow functions
def run_qiime_data_preparation(input_fp, output_dir, command_handler,\
    params, qiime_config, parallel=False,\
    status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime 
    
        The steps performed by this function are:
          1) Pick OTUs with cdhit;
          2) Pick a representative set;
          3) Align the representative set; 
          4) Assign taxonomy;
          5) Filter the alignment prior to tree building - remove positions
             which are all gaps, and specified as 0 in the lanemask
          6) Build a phylogenetic tree;
          7) Build an OTU table.
    
    """
    
    # Prepare some variables for the later steps
    input_dir, input_filename = split(input_fp)
    input_basename, input_ext = splitext(input_filename)
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    qiime_home = qiime_config['qiime_home']
    qiime_dir = qiime_config['qiime_dir']
    
    # Prep the OTU picking command
    pick_otu_dir = '%s/%s_picked_otus' % \
     (output_dir, params['pick_otus']['otu_picking_method'])
    otu_fp = '%s/%s_otus.txt' % (pick_otu_dir,input_basename)
    try:
        params_str = get_params_str(params['pick_otus'])
    except KeyError:
        params_str = ''
    # Build the OTU picking command
    pick_otus_cmd = '%s %s/pick_otus.py -i %s -o %s %s' %\
     (python_exe_fp, qiime_dir, input_fp, pick_otu_dir, params_str)
    commands.append([('Pick OTUs', pick_otus_cmd)])
    
    # Prep the representative set picking command
    rep_set_dir = '%s/rep_set/' % pick_otu_dir
    try:
        makedirs(rep_set_dir)
    except OSError:
        pass
    rep_set_fp = '%s/%s_rep_set.fasta' % (rep_set_dir,input_basename)
    rep_set_log_fp = '%s/%s_pick_rep_set.log' % (rep_set_dir,input_basename)
    try:
        params_str = get_params_str(params['pick_rep_set'])
    except KeyError:
        params_str = ''
    # Build the representative set picking command
    pick_rep_set_cmd = '%s %s/pick_rep_set.py -i %s -f %s -l %s -o %s %s' %\
     (python_exe_fp, qiime_dir, otu_fp, input_fp, rep_set_log_fp,\
      rep_set_fp, params_str)
    commands.append([('Pick representative set', pick_rep_set_cmd)])
    
    # # Set script file paths based on whether the run is in parallel or not
    # if parallel:
    #     align_seqs_fp = '%s/parallel/align_seqs_pynast.py' % qiime_dir
    #     if rdp:
    #         assign_taxonomy_fp = '%s/parallel/assign_taxonomy_rdp.py' \
    #          % qiime_dir
    #     else:
    #         assign_taxonomy_fp = '%s/parallel/assign_taxonomy_blast.py' \
    #          % qiime_dir
    # else:
    align_seqs_fp = '%s/align_seqs.py -m pynast' % qiime_dir
    assign_taxonomy_fp = '%s/assign_taxonomy.py' % qiime_dir
    
    # Prep the pynast alignment command
    pynast_dir = '%s/%s_aligned_seqs' % \
     (rep_set_dir,params['align_seqs']['alignment_method'])
    aln_fp = '%s/%s_rep_set_aligned.fasta' % (pynast_dir,input_basename)
    try:
        params_str = get_params_str(params['align_seqs'])
    except KeyError:
        params_str = ''
    # Build the pynast alignment command
    align_seqs_cmd = '%s %s -i %s -o %s %s' %\
     (python_exe_fp, align_seqs_fp, rep_set_fp, pynast_dir, params_str)
    
    # Prep the taxonomy assignment command
    assign_taxonomy_dir = '%s/%s_assigned_taxonomy' %\
     (rep_set_dir,params['assign_taxonomy']['assignment_method'])
    taxonomy_fp = '%s/%s_rep_set_tax_assignments.txt' % \
     (assign_taxonomy_dir,input_basename)
    try:
        params_str = get_params_str(params['assign_taxonomy'])
    except KeyError:
        params_str = ''
    # Build the taxonomy assignment command
    assign_taxonomy_cmd = '%s %s -o %s -i %s %s' %\
     (python_exe_fp, assign_taxonomy_fp, assign_taxonomy_dir,\
      rep_set_fp, params_str)
    
    # Append commands which can be run simulataneously in parallel
    commands.append([('Align sequences', align_seqs_cmd),\
                     ('Assign taxonomy',assign_taxonomy_cmd)])
    
    # Prep the alignment filtering command
    filtered_aln_fp = '%s/%s_rep_set_aligned_pfiltered.fasta' %\
     (pynast_dir,input_basename)
    try:
        params_str = get_params_str(params['filter_alignment'])
    except KeyError:
        params_str = ''
    # Build the alignment filtering command
    filter_alignment_cmd = '%s %s/filter_alignment.py -o %s -i %s %s' %\
     (python_exe_fp, qiime_dir, pynast_dir, aln_fp, params_str)
    commands.append([('Filter alignment', filter_alignment_cmd)])
    
    # Prep the tree building command
    phylogeny_dir = '%s/%s_phylogeny' %\
     (pynast_dir, params['make_phylogeny']['tree_method'])
    try:
        makedirs(phylogeny_dir)
    except OSError:
        pass
    tree_fp = '%s/%s_rep_set.tre' % (phylogeny_dir,input_basename)
    log_fp = '%s/%s_rep_set_phylogeny.log' % (phylogeny_dir,input_basename)
    try:
        params_str = get_params_str(params['make_phylogeny'])
    except KeyError:
        params_str = ''
    # Build the tree building command
    make_phylogeny_cmd = '%s %s/make_phylogeny.py -i %s -o %s -l %s %s' %\
     (python_exe_fp, qiime_dir, filtered_aln_fp, tree_fp, log_fp,\
     params_str)
    
    # Prep the OTU table building command
    otu_table_dir = '%s/otu_table/' % assign_taxonomy_dir
    try:
        makedirs(otu_table_dir)
    except OSError:
        pass
    otu_table_fp = '%s/%s_otu_table.txt' % (otu_table_dir,input_basename)
    try:
        params_str = get_params_str(params['make_otu_table'])
    except KeyError:
        params_str = ''
    # Build the OTU table building command
    make_otu_table_cmd = '%s %s/make_otu_table.py -i %s -t %s -o %s %s' %\
     (python_exe_fp, qiime_dir, otu_fp, taxonomy_fp, otu_table_fp, params_str)
    
    
    # Append commands which can be run simulataneously in parallel
    commands.append([('Build phylogenetic tree', make_phylogeny_cmd),\
                     ('Make OTU table', make_otu_table_cmd)])
    
    # Call the command handler on the list of commands
    command_handler(commands,status_update_callback)
## End task-specific workflow functions
    