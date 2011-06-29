#!/usr/bin/env python
# File created on 30 Dec 2009.
from __future__ import division
import sys
import re
from subprocess import Popen, PIPE, STDOUT
from os import makedirs, listdir
from glob import glob
from os.path import split, splitext, join, dirname, abspath
from datetime import datetime
from numpy import array
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import parse_mapping_file, parse_qiime_parameters
from qiime.format import format_otu_table
from qiime.util import (compute_seqs_per_library_stats,
                        get_qiime_scripts_dir,
                        create_dir, guess_even_sampling_depth,
                        get_interesting_mapping_fields,qiime_system_call,
                        get_qiime_library_version)
from cogent.core.moltype import IUPAC_DNA_ambiguities
import os

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

"""
This file contains the QIIME workflow functions which string together 
independent scripts. For usage examples see the related files in the 
scripts directory:
 - 
"""

## Start utilities used by the workflow functions
def generate_log_fp(output_dir,
                    basefile_name='log',
                    suffix='txt',
                    timestamp_pattern='%Y%m%d%H%M%S'):
    timestamp = datetime.now().strftime(timestamp_pattern)
    filename = '%s_%s.%s' % (basefile_name,timestamp,suffix)
    return join(output_dir,filename)

class WorkflowError(Exception):
    pass

class WorkflowLogger(object):
    
    def __init__(self,log_fp=None,params=None,qiime_config=None,open_mode='w'):
        if log_fp:
            self._f = open(log_fp,open_mode)
        else:
            self._f = None
        start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('Logging started at %s\n' % start_time)
        self.write('QIIME version: %s\n\n' % get_qiime_library_version())
        self.writeQiimeConfig(qiime_config)
        self.writeParams(params)
    
    def write(self,s):
        if self._f:
            self._f.write(s)
            # Flush here so users can see what step they're
            # on after each write, since some steps can take
            # a long time, and a relatively small amount of 
            # data is being written to the log files.
            self._f.flush()
        else:
            pass
    
    def writeQiimeConfig(self,qiime_config):
        if qiime_config == None:
            self.write('No qiime config provided.\n')
        else:
            self.write('qiime_config values:\n')
            for k,v in qiime_config.items():
                if v:
                    self.write('%s\t%s\n' % (k,v))
            self.write('\n')
            
    def writeParams(self,params):
        if params == None:
            self.write('No params provided.\n')
        else:
            self.write('parameter file values:\n')
            for k,v in params.items():
                for inner_k,inner_v in v.items():
                    val = inner_v or 'True'
                    self.write('%s:%s\t%s\n' % (k,inner_k,val))
            self.write('\n')
    
    def close(self):
        end_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('\nLogging stopped at %s\n' % end_time)
        if self._f:
            self._f.close()
        else:
            pass

def print_commands(commands,
                   status_update_callback,
                   logger,
                   close_logger_on_success=True):
    """Print list of commands to run """
    logger.write("Printing commands only.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('#%s' % e[0])
            print '%s' % e[1]
            logger.write('# %s command\n%s\n\n' % e)
            
def call_commands_serially(commands,
                           status_update_callback,
                           logger,
                           close_logger_on_success=True):
    """Run list of commands, one after another """
    logger.write("Executing commands.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('%s\n%s' % e)
            logger.write('# %s command \n%s\n\n' % e)
            stdout, stderr, return_value = qiime_system_call(e[1])
            if return_value != 0:
                msg = "\n\n*** ERROR RAISED DURING STEP: %s\n" % e[0] +\
                 "Command run was:\n %s\n" % e[1] +\
                 "Command returned exit status: %d\n" % return_value +\
                 "Stdout:\n%s\nStderr\n%s\n" % (stdout,stderr)
                logger.write(msg)
                logger.close()
                raise WorkflowError, msg
            # in the no error case, we write commands' output to the log
            # and also echo to this proc's stdout/stderr
            else:
                # write stdout and stderr to log file
                logger.write("Stdout:\n%s\nStderr:\n%s\n" % (stdout,stderr))
                # write stdout to stdout
                if stdout:
                    print stdout
                # write stderr to stderr
                if stderr:
                    sys.stderr.write(stderr)
    if close_logger_on_success: logger.close()

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

def validate_and_set_jobs_to_start(params,
                                   jobs_to_start,
                                   default_jobs_to_start,
                                   parallel,
                                   option_parser):
    if (jobs_to_start != int(default_jobs_to_start)) and \
       not parallel:
        option_parser.error("Passing -O requires that -a is also passed.")
    params['parallel']['jobs_to_start'] = str(jobs_to_start)
    

## End utilities used by the workflow functions

## Begin task-specific workflow functions
def run_qiime_data_preparation(input_fp, 
                               output_dir, 
                               command_handler,
                               params, 
                               qiime_config,
                               parallel=False,
                               logger=None,
                               status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime 
    
        The steps performed by this function are:
          1) Pick OTUs;
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
    create_dir(output_dir)
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    if logger == None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False
    
    # Prep the OTU picking command
    try:
        otu_picking_method = params['pick_otus']['otu_picking_method']
    except KeyError:
        otu_picking_method = 'uclust'
    pick_otu_dir = '%s/%s_picked_otus' % (output_dir, otu_picking_method)
    otu_fp = '%s/%s_otus.txt' % (pick_otu_dir,input_basename)
    if parallel and (otu_picking_method == 'blast' or 
                     otu_picking_method == 'uclust_ref'):
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''
        
        # Grab the OTU picker parameters
        try:
            # Want to find a cleaner strategy for this: the parallel script
            # is method-specific, so doesn't take a --otu_picking_method
            # option. This works for now though.
            d = params['pick_otus'].copy()
            del d['otu_picking_method']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass
        otu_picking_script = 'parallel_pick_otus_%s.py' % otu_picking_method
        # Build the OTU picking command
        pick_otus_cmd = '%s %s/%s -i %s -o %s -T %s' % (python_exe_fp, 
                                                        script_dir, 
                                                        otu_picking_script,
                                                        input_fp,
                                                        pick_otu_dir,
                                                        params_str)
    else:
        try:
            params_str = get_params_str(params['pick_otus'])
        except KeyError:
            params_str = ''
        # Build the OTU picking command
        pick_otus_cmd = '%s %s/pick_otus.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, input_fp, pick_otu_dir, params_str)

    commands.append([('Pick OTUs', pick_otus_cmd)])
    
    # Prep the representative set picking command
    rep_set_dir = '%s/rep_set/' % output_dir
    try:
        makedirs(rep_set_dir)
    except OSError:
        pass
    rep_set_fp = '%s/%s_rep_set.fasta' % (rep_set_dir,input_basename)
    rep_set_log_fp = '%s/%s_rep_set.log' % (rep_set_dir,input_basename)
    
    try:
        params_str = get_params_str(params['pick_rep_set'])
    except KeyError:
        params_str = ''
    # Build the representative set picking command
    pick_rep_set_cmd = '%s %s/pick_rep_set.py -i %s -f %s -l %s -o %s %s' %\
     (python_exe_fp, script_dir, otu_fp, input_fp, rep_set_log_fp,\
      rep_set_fp, params_str)
    commands.append([('Pick representative set', pick_rep_set_cmd)])
    
    # Prep the taxonomy assignment command
    try:
        assignment_method = params['assign_taxonomy']['assignment_method']
    except KeyError:
        assignment_method = 'rdp'
    assign_taxonomy_dir = '%s/%s_assigned_taxonomy' %\
     (output_dir,assignment_method)
    taxonomy_fp = '%s/%s_rep_set_tax_assignments.txt' % \
     (assign_taxonomy_dir,input_basename)
    if parallel and (assignment_method == 'rdp' or assignment_method == 'blast'):
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''
        
        # Grab the OTU picker parameters
        try:
            # Want to find a cleaner strategy for this: the parallel script
            # is method-specific, so doesn't take a --assignment_method
            # option. This works for now though.
            d = params['assign_taxonomy'].copy()
            del d['assignment_method']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass
            
        # Build the parallel taxonomy assignment command
        assign_taxonomy_cmd = \
         '%s %s/parallel_assign_taxonomy_%s.py -i %s -o %s -T %s' %\
         (python_exe_fp, script_dir, assignment_method, rep_set_fp,\
          assign_taxonomy_dir, params_str)
    else:
        try:
            params_str = get_params_str(params['assign_taxonomy'])
        except KeyError:
            params_str = ''
        # Build the taxonomy assignment command
        assign_taxonomy_cmd = '%s %s/assign_taxonomy.py -o %s -i %s %s' %\
         (python_exe_fp, script_dir, assign_taxonomy_dir,\
          rep_set_fp, params_str)
    
    commands.append([('Assign taxonomy',assign_taxonomy_cmd)])
    
    # Prep the OTU table building command
    otu_table_fp = '%s/otu_table.txt' % output_dir
    try:
        params_str = get_params_str(params['make_otu_table'])
    except KeyError:
        params_str = ''
    # Build the OTU table building command
    make_otu_table_cmd = '%s %s/make_otu_table.py -i %s -t %s -o %s %s' %\
     (python_exe_fp, script_dir, otu_fp, taxonomy_fp, otu_table_fp, params_str)
    
    commands.append([('Make OTU table', make_otu_table_cmd)])
    
    # Prep the pynast alignment command
    try:
        alignment_method = params['align_seqs']['alignment_method']
    except KeyError:
        alignment_method = 'pynast'
    pynast_dir = '%s/%s_aligned_seqs' % (output_dir,alignment_method)
    aln_fp = '%s/%s_rep_set_aligned.fasta' % (pynast_dir,input_basename)
    if parallel and alignment_method == 'pynast':
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''
        
        # Grab the OTU picker parameters
        try:
            # Want to find a cleaner strategy for this: the parallel script
            # is method-specific, so doesn't take a --alignment_method
            # option. This works for now though.
            d = params['align_seqs'].copy()
            del d['alignment_method']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass
            
        # Build the parallel pynast alignment command
        align_seqs_cmd = '%s %s/parallel_align_seqs_pynast.py -i %s -o %s -T %s' %\
         (python_exe_fp, script_dir, rep_set_fp, pynast_dir, params_str)
    else:
        try:
            params_str = get_params_str(params['align_seqs'])
        except KeyError:
            params_str = ''
        # Build the pynast alignment command
        align_seqs_cmd = '%s %s/align_seqs.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, rep_set_fp, pynast_dir, params_str)
    commands.append([('Align sequences', align_seqs_cmd)])
    
    if alignment_method == 'pynast':
        # Prep the alignment filtering command (only applicable when aligned
        # with pynast)
        filtered_aln_fp = '%s/%s_rep_set_aligned_pfiltered.fasta' %\
         (pynast_dir,input_basename)
        try:
            params_str = get_params_str(params['filter_alignment'])
        except KeyError:
            params_str = ''
        # Build the alignment filtering command
        filter_alignment_cmd = '%s %s/filter_alignment.py -o %s -i %s %s' %\
         (python_exe_fp, script_dir, pynast_dir, aln_fp, params_str)
        commands.append([('Filter alignment', filter_alignment_cmd)])
    else: 
        filtered_aln_fp = aln_fp
    
    # Prep the tree building command
    tree_fp = '%s/rep_set.tre' % output_dir
    try:
        params_str = get_params_str(params['make_phylogeny'])
    except KeyError:
        params_str = ''
    # Build the tree building command
    make_phylogeny_cmd = '%s %s/make_phylogeny.py -i %s -o %s %s' %\
     (python_exe_fp, script_dir, filtered_aln_fp, tree_fp,\
     params_str)
    commands.append([('Build phylogenetic tree', make_phylogeny_cmd)])
    
    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)
    
    return abspath(tree_fp), abspath(otu_table_fp)
    
    
    
## Start reference otu picking workflow


## Begin task-specific workflow functions
def run_pick_reference_otus_through_otu_table(
                              input_fp, 
                              refseqs_fp,
                              output_dir,
                              taxonomy_fp,
                              command_handler,
                              params,
                              qiime_config,
                              parallel=False,
                              logger=None,
                              status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime 
    
        The steps performed by this function are:
          1) Pick OTUs;
          2) Build an OTU table with optional pre-defined taxonmy.
    
    """
    
    # confirm that a valid otu picking method was supplied before doing
    # any work
    reference_otu_picking_methods = ['blast','uclust_ref']

    try:
        otu_picking_method = params['pick_otus']['otu_picking_method']
    except KeyError:
        otu_picking_method = 'uclust_ref'
    assert otu_picking_method in reference_otu_picking_methods,\
     "Invalid OTU picking method supplied: %s. Valid choices are: %s"\
     % (otu_picking_method,' '.join(reference_otu_picking_methods))
    
    # Prepare some variables for the later steps
    input_dir, input_filename = split(input_fp)
    input_basename, input_ext = splitext(input_filename)
    create_dir(output_dir)
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    if logger == None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False
    
    # Prep the OTU picking command
    pick_otu_dir = '%s/%s_picked_otus' % (output_dir, otu_picking_method)
    otu_fp = '%s/%s_otus.txt' % (pick_otu_dir,input_basename)
    if parallel and (otu_picking_method == 'blast' or 
                     otu_picking_method == 'uclust_ref'):
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''
        
        # Grab the OTU picker parameters
        try:
            # Want to find a cleaner strategy for this: the parallel script
            # is method-specific, so doesn't take a --otu_picking_method
            # option. This works for now though.
            d = params['pick_otus'].copy()
            del d['otu_picking_method']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass
        otu_picking_script = 'parallel_pick_otus_%s.py' % otu_picking_method
        # Build the OTU picking command
        pick_otus_cmd = '%s %s/%s -i %s -o %s -r %s -T %s' %\
          (python_exe_fp, 
           script_dir, 
           otu_picking_script,
           input_fp,
           pick_otu_dir,
           refseqs_fp,
           params_str)
    else:
        try:
            params_str = get_params_str(params['pick_otus'])
        except KeyError:
            params_str = ''
        # Build the OTU picking command
        pick_otus_cmd = '%s %s/pick_otus.py -i %s -o %s -r %s %s' %\
         (python_exe_fp,
          script_dir,
          input_fp,
          pick_otu_dir,
          refseqs_fp,
          params_str)

    commands.append([('Pick OTUs', pick_otus_cmd)])

    # Prep the OTU table building command
    otu_table_fp = '%s/otu_table.txt' % pick_otu_dir
    try:
        params_str = get_params_str(params['make_otu_table'])
    except KeyError:
        params_str = ''
    if taxonomy_fp:
        taxonomy_str = '-t %s' % taxonomy_fp
    else:
        taxonomy_str = ''
    # Build the OTU table building command
    make_otu_table_cmd = '%s %s/make_otu_table.py -i %s %s -o %s %s' %\
     (python_exe_fp, script_dir, otu_fp, taxonomy_str, otu_table_fp, params_str)
    
    commands.append([('Make OTU table', make_otu_table_cmd)])
    

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)



def run_beta_diversity_through_plots(otu_table_fp, mapping_fp,
    output_dir, command_handler, params, qiime_config,
    color_by_interesting_fields_only=True,sampling_depth=None,
    tree_fp=None, parallel=False, logger=None, suppress_3d_plots=False,
    suppress_2d_plots=False, suppress_distance_histograms=False,
    status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime 
    
        The steps performed by this function are:
         1) Compute a beta diversity distance matrix;
         2) Peform a principal coordinates analysis on the result of
          Step 1;
         3) Generate a 3D prefs file for optimized coloring of continuous
          variables;
         4) Generate a 3D plot for all mapping fields with colors
          optimized for continuous data;
         5) Generate a 3D plot for all mapping fields with colors
          optimized for discrete data.
    
    """  
    # Prepare some variables for the later steps
    otu_table_dir, otu_table_filename = split(otu_table_fp)
    otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    create_dir(output_dir)
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    if logger == None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False
    
    mapping_data, mapping_header, mapping_comments =\
     parse_mapping_file(open(mapping_fp,'U'))
    # Get the interesting mapping fields to color by -- if none are
    # interesting, take all of them. Interesting is defined as those
    # which have greater than one value and fewer values than the number 
    # of samples
    if color_by_interesting_fields_only:
        mapping_fields =\
          get_interesting_mapping_fields(mapping_data, mapping_header) or mapping_header
    else:
        mapping_fields = mapping_header
    mapping_fields = ','.join(mapping_fields)
    
    if sampling_depth:
        # Sample the OTU table at even depth
        even_sampled_otu_table_fp = '%s/%s_even%d%s' %\
         (output_dir, otu_table_basename, 
          sampling_depth, otu_table_ext)
        single_rarefaction_cmd = \
         '%s %s/single_rarefaction.py -i %s -o %s -d %d' %\
         (python_exe_fp, script_dir, otu_table_fp,
          even_sampled_otu_table_fp, sampling_depth)
        commands.append([
         ('Sample OTU table at %d seqs/sample' % sampling_depth,
          single_rarefaction_cmd)])
        otu_table_fp = even_sampled_otu_table_fp
        otu_table_dir, otu_table_filename = split(even_sampled_otu_table_fp)
        otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    try:
        beta_diversity_metrics = params['beta_diversity']['metrics'].split(',')
    except KeyError:
        beta_diversity_metrics = ['weighted_unifrac','unweighted_unifrac']
    
    # Prep the 3d prefs file generator command
    prefs_fp = '%s/prefs.txt' % output_dir
    try:
        params_str = get_params_str(params['make_prefs_file'])
    except KeyError:
        params_str = ''
    if not 'mapping_headers_to_use' in params['make_prefs_file']:
        params_str = '%s --mapping_headers_to_use %s' \
         % (params_str,mapping_fields)
    # Build the 3d prefs file generator command
    prefs_cmd = \
     '%s %s/make_prefs_file.py -m %s -o %s %s' %\
     (python_exe_fp, script_dir, mapping_fp, prefs_fp, params_str)
    commands.append([('Build prefs file', prefs_cmd)])
    
    dm_fps = []
    for beta_diversity_metric in beta_diversity_metrics:
        
        # Prep the beta-diversity command
        try:
            bdiv_params_copy = params['beta_diversity'].copy()
        except KeyError:
            bdiv_params_copy = {}
        try:
            del bdiv_params_copy['metrics']
        except KeyError:
            pass
        
        params_str = get_params_str(bdiv_params_copy)
            
        if tree_fp:
            params_str = '%s -t %s ' % (params_str,tree_fp)
            
        # Build the beta-diversity command
        if parallel:
            # Grab the parallel-specific parameters
            try:
                params_str += get_params_str(params['parallel'])
            except KeyError:
                pass
            beta_div_cmd = '%s %s/parallel_beta_diversity.py -i %s -o %s --metrics %s -T %s' %\
             (python_exe_fp, script_dir, otu_table_fp,
              output_dir, beta_diversity_metric, params_str)
            commands.append(\
             [('Beta Diversity (%s)' % beta_diversity_metric, beta_div_cmd)])
        else:
            beta_div_cmd = '%s %s/beta_diversity.py -i %s -o %s --metrics %s %s' %\
             (python_exe_fp, script_dir, otu_table_fp, 
              output_dir, beta_diversity_metric, params_str)
            commands.append(\
             [('Beta Diversity (%s)' % beta_diversity_metric, beta_div_cmd)])
        
        
        orig_beta_div_fp = '%s/%s_%s' % \
         (output_dir, beta_diversity_metric, otu_table_filename)
        beta_div_fp = '%s/%s_dm.txt' % \
         (output_dir, beta_diversity_metric)
        commands.append([('Rename distance matrix (%s)' % beta_diversity_metric,
                         'mv %s %s' % (orig_beta_div_fp, beta_div_fp))])
        dm_fps.append((beta_diversity_metric, beta_div_fp))
        
        # Prep the principal coordinates command
        pc_fp = '%s/%s_pc.txt' % (output_dir, beta_diversity_metric)
        try:
            params_str = get_params_str(params['principal_coordinates'])
        except KeyError:
            params_str = ''
        # Build the principal coordinates command
        pc_cmd = '%s %s/principal_coordinates.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, beta_div_fp, pc_fp, params_str)
        commands.append(\
         [('Principal coordinates (%s)' % beta_diversity_metric, pc_cmd)])
        
        # Generate 3d plots
        if not suppress_3d_plots:
            # Prep the continuous-coloring 3d plots command
            continuous_3d_dir = '%s/%s_3d_continuous/' %\
             (output_dir, beta_diversity_metric)
            try:
                makedirs(continuous_3d_dir)
            except OSError:
                pass
            try:
                params_str = get_params_str(params['make_3d_plots'])
            except KeyError:
                params_str = ''
            # Build the continuous-coloring 3d plots command
            continuous_3d_command = \
             '%s %s/make_3d_plots.py -p %s -i %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, prefs_fp, pc_fp, continuous_3d_dir,\
               mapping_fp, params_str)
    
            # Prep the discrete-coloring 3d plots command
            discrete_3d_dir = '%s/%s_3d_discrete/' %\
             (output_dir, beta_diversity_metric)
            try:
                makedirs(discrete_3d_dir)
            except OSError:
                pass
            try:
                params_str = get_params_str(params['make_3d_plots'])
            except KeyError:
                params_str = ''
            # Build the discrete-coloring 3d plots command
            discrete_3d_command = \
             '%s %s/make_3d_plots.py -b "%s" -i %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, mapping_fields, pc_fp, discrete_3d_dir,\
               mapping_fp, params_str)
       
            commands.append([\
              ('Make 3D plots (continuous coloring, %s)' %\
                beta_diversity_metric,continuous_3d_command),\
              ('Make 3D plots (discrete coloring, %s)' %\
                beta_diversity_metric,discrete_3d_command,)])
    
        # Generate 3d plots
        if not suppress_2d_plots:
            # Prep the continuous-coloring 3d plots command
            continuous_2d_dir = '%s/%s_2d_continuous/' %\
             (output_dir, beta_diversity_metric)
            try:
                makedirs(continuous_2d_dir)
            except OSError:
                pass
            try:
                params_str = get_params_str(params['make_2d_plots'])
            except KeyError:
                params_str = ''
            # Build the continuous-coloring 3d plots command
            continuous_2d_command = \
             '%s %s/make_2d_plots.py -p %s -i %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, prefs_fp, pc_fp, continuous_2d_dir,\
               mapping_fp, params_str)
               
            # Prep the discrete-coloring 3d plots command
            discrete_2d_dir = '%s/%s_2d_discrete/' %\
             (output_dir, beta_diversity_metric)
            try:
                makedirs(discrete_2d_dir)
            except OSError:
                pass
            try:
                params_str = get_params_str(params['make_2d_plots'])
            except KeyError:
                params_str = ''
            # Build the discrete-coloring 2d plots command
            discrete_2d_command = \
             '%s %s/make_2d_plots.py -b "%s" -i %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, mapping_fields, pc_fp, discrete_2d_dir,\
               mapping_fp, params_str)
       
            commands.append([\
              ('Make 2D plots (continuous coloring, %s)' %\
                beta_diversity_metric,continuous_2d_command),\
              ('Make 2D plots (discrete coloring, %s)' %\
                beta_diversity_metric,discrete_2d_command,)])
                
        if not suppress_distance_histograms:
            # Prep the discrete-coloring 3d plots command
            histograms_dir = '%s/%s_histograms/' %\
             (output_dir, beta_diversity_metric)
            try:
                makedirs(histograms_dir)
            except OSError:
                pass
            try:
                params_str = get_params_str(params['make_distance_histograms'])
            except KeyError:
                params_str = ''
            # Build the make_distance_histograms command
            distance_histograms_command = \
             '%s %s/make_distance_histograms.py -d %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, beta_div_fp, 
               histograms_dir, mapping_fp, params_str)
       
            commands.append([\
              ('Make Distance Histograms (%s)' %\
                beta_diversity_metric,distance_histograms_command)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)
    
    return dm_fps


def run_qiime_alpha_rarefaction(otu_table_fp, mapping_fp,
    output_dir, command_handler, params, qiime_config, tree_fp=None,
    num_steps=10, parallel=False, logger=None, min_seqs_per_sample=10,
    status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime 
    
        The steps performed by this function are:
          1) Generate rarefied OTU tables;
          2) Compute alpha diversity metrics for each rarefied OTU table;
          3) Collate alpha diversity results;
          4) Generate alpha rarefaction plots.
    
    """
    # Prepare some variables for the later steps
    otu_table_dir, otu_table_filename = split(otu_table_fp)
    otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    create_dir(output_dir)
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    if logger == None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False
    
    # Prep the rarefaction command
    try:
        otu_table_f = open(otu_table_fp,'U')
    except IOError,e:
        logger.write('OTU table filepath cannot be opened. Does it exist?\n' +
                     ' %s\n' % otu_table_fp +
                     'Original Error:\n%s\n' % str(e))
        logger.close()
        raise IOError,e
    
    min_count, max_count, median_count, mean_count, counts_per_sample =\
     compute_seqs_per_library_stats(otu_table_f)
    step = int((median_count - min_seqs_per_sample) / num_steps) or 1
    median_count = int(median_count)
    
    rarefaction_dir = '%s/rarefaction/' % output_dir
    try:
        makedirs(rarefaction_dir)
    except OSError:
        pass
    try:
        params_str = get_params_str(params['multiple_rarefactions'])
    except KeyError:
        params_str = ''
    if parallel:
        params_str += ' %s' % get_params_str(params['parallel'])        
        # Build the rarefaction command
        rarefaction_cmd = \
         '%s %s/parallel_multiple_rarefactions.py -T -i %s -m %s -x %s -s %s -o %s %s' %\
         (python_exe_fp, script_dir, otu_table_fp, min_seqs_per_sample, median_count, \
          step, rarefaction_dir, params_str)
    else:
        # Build the rarefaction command
        rarefaction_cmd = \
         '%s %s/multiple_rarefactions.py -i %s -m %s -x %s -s %s -o %s %s' %\
         (python_exe_fp, script_dir, otu_table_fp, min_seqs_per_sample, median_count, \
          step, rarefaction_dir, params_str)
    commands.append([('Alpha rarefaction', rarefaction_cmd)])
    
    # Prep the alpha diversity command
    alpha_diversity_dir = '%s/alpha_div/' % output_dir
    try:
        makedirs(alpha_diversity_dir)
    except OSError:
        pass
    try:
        params_str = get_params_str(params['alpha_diversity'])
    except KeyError:
        params_str = ''
    if parallel:
        params_str += ' %s' % get_params_str(params['parallel'])   
        # Build the alpha diversity command
        alpha_diversity_cmd = \
         "%s %s/parallel_alpha_diversity.py -T -i %s -o %s -t %s %s" %\
         (python_exe_fp, script_dir, rarefaction_dir, alpha_diversity_dir, \
          tree_fp, params_str)
    else:  
        # Build the alpha diversity command
        alpha_diversity_cmd = \
         "%s %s/alpha_diversity.py -i %s -o %s -t %s %s" %\
         (python_exe_fp, script_dir, rarefaction_dir, alpha_diversity_dir, \
          tree_fp, params_str)

    commands.append(\
     [('Alpha diversity on rarefied OTU tables',alpha_diversity_cmd)])
     
    # Prep the alpha diversity collation command
    # python $qdir/collate_alpha.py -i Fasting_Alpha_Metrics/ -o Fasting_Alpha_Collated/
    alpha_collated_dir = '%s/alpha_div_collated/' % output_dir
    try:
        makedirs(alpha_collated_dir)
    except OSError:
        pass
    try:
        params_str = get_params_str(params['collate_alpha'])
    except KeyError:
        params_str = ''
    # Build the alpha diversity collation command
    alpha_collated_cmd = '%s %s/collate_alpha.py -i %s -o %s %s' %\
     (python_exe_fp, script_dir, alpha_diversity_dir, \
      alpha_collated_dir, params_str)
    commands.append([('Collate alpha',alpha_collated_cmd)])

    # Prep the make rarefaction plot command(s)
    rarefaction_plot_dir = '%s/alpha_rarefaction_plots/' % output_dir
    try:
        makedirs(rarefaction_plot_dir)
    except OSError:
        pass
    try:
        params_str = get_params_str(params['make_rarefaction_plots'])
    except KeyError:
        params_str = ''
    # Build the make rarefaction plot command(s)
    #for metric in alpha_diversity_metrics:
    make_rarefaction_plot_cmd =\
         '%s %s/make_rarefaction_plots.py -i %s -m %s -o %s %s' %\
         (python_exe_fp, script_dir, alpha_collated_dir, mapping_fp,
          rarefaction_plot_dir, params_str)
    commands.append(\
         [('Rarefaction plot: %s' % 'All metrics',make_rarefaction_plot_cmd)])
    
    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)

def run_jackknifed_beta_diversity(otu_table_fp,tree_fp,seqs_per_sample,
    output_dir, command_handler, params, qiime_config, mapping_fp,
    parallel=False,logger=None,
    status_update_callback=print_to_stdout, master_tree=None):
    """ Run the data preparation steps of Qiime 
    
        The steps performed by this function are:
          1) Compute beta diversity distance matrix from otu table (and
           tree, if applicable)
          2) Build rarefied OTU tables;
          3) Build UPGMA tree from full distance matrix;
          4) Compute distance matrics for rarefied OTU tables;
          5) Build UPGMA trees from rarefied OTU table distance matrices;
          5.5) Build a consensus tree from the rarefied UPGMA trees
          6) Compare rarefied OTU table distance matrix UPGMA trees 
           to tree full UPGMA tree and write support file and newick tree
           with support values as node labels.
           
        master_tree can be 'full' or 'consensus', default full
    """
    # Prepare some variables for the later steps
    if master_tree == None:
        master_tree = 'full'
    otu_table_dir, otu_table_filename = split(otu_table_fp)
    otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    create_dir(output_dir)
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    if logger == None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False
    try:
        beta_diversity_metrics = params['beta_diversity']['metrics'].split(',')
    except KeyError:
        beta_diversity_metrics = ['weighted_unifrac','unweighted_unifrac']
    
    # Prep the beta-diversity command
    try:
        params_str = get_params_str(params['beta_diversity'])
    except KeyError:
        params_str = ''
    if tree_fp:
        params_str = '%s -t %s' % (params_str,tree_fp)
    # Build the beta-diversity command
    beta_div_cmd = '%s %s/beta_diversity.py -i %s -o %s %s' %\
     (python_exe_fp, script_dir, otu_table_fp, output_dir, params_str)
    commands.append(\
     [('Beta Diversity (%s)' % ', '.join(beta_diversity_metrics), beta_div_cmd)])

    # Prep rarefaction command
    rarefaction_dir = '%s/rarefaction/' % output_dir
    try:
        makedirs(rarefaction_dir)
    except OSError:
        pass
    try:
        params_str = get_params_str(params['multiple_rarefactions_even_depth'])
    except KeyError:
        params_str = ''
    # if parallel:
    #     params_str += ' %s' % get_params_str(params['parallel'])  
    #     # Build the parallel rarefaction command
    #     rarefaction_cmd = \
    #      '%s %s/parallel_multiple_rarefactions.py -T -i %s -m %s -x %s -s 1 -o %s %s' %\
    #      (python_exe_fp, script_dir, otu_table_fp, seqs_per_sample,\
    #       seqs_per_sample, rarefaction_dir, params_str)
    # else:
    # Build the serial rarefaction command
    rarefaction_cmd = \
     '%s %s/multiple_rarefactions_even_depth.py -i %s -d %d -o %s %s' %\
     (python_exe_fp, script_dir, otu_table_fp, seqs_per_sample, \
      rarefaction_dir, params_str)
    commands.append([('Rarefaction', rarefaction_cmd)])

    # Begin iterating over beta diversity distance metrics, if more than one
    # was provided
    for beta_diversity_metric in beta_diversity_metrics:
        metric_output_dir = '%s/%s/' % (output_dir, beta_diversity_metric)
        distance_matrix_fp = '%s/%s_%s.txt' % \
         (output_dir, beta_diversity_metric, otu_table_basename)
    
        # Prep the hierarchical clustering command (for full distance matrix)
        full_tree_fp = '%s/%s_upgma.tre' % (metric_output_dir,otu_table_basename)
        try:
            params_str = get_params_str(params['upgma_cluster'])
        except KeyError:
            params_str = ''
        # Build the hierarchical clustering command (for full distance matrix)
        hierarchical_cluster_cmd = '%s %s/upgma_cluster.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, distance_matrix_fp, full_tree_fp, params_str)
        commands.append(\
         [('UPGMA on full distance matrix: %s' % beta_diversity_metric,\
           hierarchical_cluster_cmd)])
           
           
           
        # Prep the beta diversity command (for rarefied OTU tables)
        dm_dir = '%s/rare_dm/' % metric_output_dir
        try:
            makedirs(dm_dir)
        except OSError:
            pass
        # the metrics parameter needs to be ignored as we need to run
        # beta_diversity one metric at a time to keep the per-metric
        # output files in separate directories
        try:
            d = params['beta_diversity'].copy()
            del d['metrics']
        except KeyError:
            params_str = {}
        params_str = get_params_str(d) + ' -m %s ' % beta_diversity_metric
        if tree_fp:
            params_str = '%s -t %s' % (params_str,tree_fp)
        if parallel:
            params_str += ' %s' % get_params_str(params['parallel'])        
            # Build the parallel beta diversity command (for rarefied OTU tables)
            beta_div_rarefied_cmd = \
             '%s %s/parallel_beta_diversity.py -T -i %s -o %s %s' %\
             (python_exe_fp, script_dir, rarefaction_dir, dm_dir, params_str)
        else:
            # Build the serial beta diversity command (for rarefied OTU tables)
            beta_div_rarefied_cmd = \
             '%s %s/beta_diversity.py -i %s -o %s %s' %\
             (python_exe_fp, script_dir, rarefaction_dir, dm_dir, params_str)
        commands.append(\
         [('Beta diversity on rarefied OTU tables (%s)' % beta_diversity_metric,\
           beta_div_rarefied_cmd)])

        # Prep the hierarchical clustering command (for rarefied 
        # distance matrices)
        upgma_dir = '%s/rare_upgma/' % metric_output_dir
        try:
            makedirs(upgma_dir)
        except OSError:
            pass

        try:
            params_str = get_params_str(params['upgma_cluster'])
        except KeyError:
            params_str = ''
        # Build the hierarchical clustering command (for rarefied 
        # distance matrices)
        hierarchical_cluster_cmd =\
         '%s %s/upgma_cluster.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, dm_dir, upgma_dir, params_str)
        commands.append(\
         [('UPGMA on rarefied distance matrix (%s)' % beta_diversity_metric,\
           hierarchical_cluster_cmd)])
        

        # Build the consensus tree command
        consensus_tree_cmd =\
         '%s %s/consensus_tree.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, upgma_dir, upgma_dir + "/consensus.tre",
            params_str)
        commands.append(\
         [('consensus on rarefied distance matrices (%s)' % beta_diversity_metric,\
           consensus_tree_cmd)])
           
           
        # Prep the tree compare command
        tree_compare_dir = '%s/upgma_cmp/' % metric_output_dir
        try:
            makedirs(tree_compare_dir)
        except OSError:
            pass
        try:
            params_str = get_params_str(params['tree_compare'])
        except KeyError:
            params_str = ''

        # Build the tree compare command
        if master_tree == "full":
            master_tree_fp = full_tree_fp
        elif master_tree == "consensus":
            master_tree_fp = upgma_dir + "/consensus.tre"
        else:
            raise RuntimeError('master tree method "%s" not found' % (master_tree,))
        tree_compare_cmd = '%s %s/tree_compare.py -s %s -m %s -o %s %s' %\
         (python_exe_fp, script_dir, upgma_dir, master_tree_fp, \
          tree_compare_dir, params_str)
        commands.append(\
         [('Tree compare (%s)' % beta_diversity_metric,\
           tree_compare_cmd)])
           
        # Prep the PCoA command
        pcoa_dir = '%s/pcoa/' % metric_output_dir
        try:
            makedirs(pcoa_dir)
        except OSError:
            pass
        try:
            params_str = get_params_str(params['principal_coordinates'])
        except KeyError:
            params_str = ''
        # Build the PCoA command
        pcoa_cmd = '%s %s/principal_coordinates.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, dm_dir, pcoa_dir, params_str)
        commands.append(\
         [('Principal coordinates (%s)' % beta_diversity_metric, pcoa_cmd)])
           
        # Prep the 2D plots command
        plots_2d_dir = '%s/2d_plots/' % metric_output_dir
        try:
            makedirs(plots_2d_dir)
        except OSError:
            pass
        try:
            params_str = get_params_str(params['make_2d_plots'])
        except KeyError:
            params_str = ''
        # Build the 2d plots command
        plots_2d_cmd = '%s %s/make_2d_plots.py -i %s -o %s -m %s %s' %\
         (python_exe_fp, script_dir, pcoa_dir, plots_2d_dir, 
          mapping_fp, params_str)
        commands.append(\
         [('2d plots (%s)' % beta_diversity_metric, plots_2d_cmd)])
         
        # Prep the 3D plots command
        plots_3d_dir = '%s/3d_plots/' % metric_output_dir
        try:
            makedirs(plots_3d_dir)
        except OSError:
            pass
        try:
            params_str = get_params_str(params['make_3d_plots'])
        except KeyError:
            params_str = ''
        # Build the 2d plots command
        plots_3d_cmd = '%s %s/make_3d_plots.py -i %s -o %s -m %s %s' %\
         (python_exe_fp, script_dir, pcoa_dir, plots_3d_dir, 
          mapping_fp, params_str)
        commands.append(\
         [('3d plots (%s)' % beta_diversity_metric, plots_3d_cmd)])
           
           

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)
    


def format_index_link(link_description,relative_path):
    
    return '<td>%s</td><td> <a href="%s">%s</a></td>' % (link_description,
                                        re.sub('/+','/',relative_path),
                                        split(relative_path)[1])

def generate_index_page(index_links,index_fp):
    # get containing directory for index_fp
    top_level_dir = split(split(index_fp)[0])[1]
    index_page_header = "<html><head><title>QIIME results</title></head><body>\n"
    index_lines = [index_page_header]
    d = {}
    for e in index_links:
        try:
            d[e[2]].append((e[0],e[1]))
        except KeyError:
            d[e[2]] = [(e[0],e[1])]
    index_lines.append('<table border=1>\n')
    for k,v in d.items():
        index_lines.append(
         '<tr colspan=2 align=center bgcolor=wheat><td colspan=2 align=center>%s</td></tr>\n' % k)
        for description,path in v:
            path = re.sub('.*%s' % top_level_dir,'./',path)
            index_lines.append('<tr>%s</tr>\n' % format_index_link(description,path))
    index_lines.append('</table>')
    
    index_page_footer = "</body></html>"
    index_lines.append(index_page_footer)
    
    open(index_fp,'w').write(''.join(index_lines))

# Run QIIME method comparison workflow
def run_core_qiime_analyses(
    fna_fps,
    qual_fps,
    mapping_fp,
    output_dir,
    command_handler,
    qiime_config,
    params=None,
    categories=None,
    sampling_depth=None,
    even_sampling_keeps_all_samples=False,
    suppress_split_libraries=False,
    arare_min_seqs_per_sample=10,
    arare_num_steps=10,
    reference_tree_fp=None,
    parallel=False,
    status_update_callback=print_to_stdout):
    """ Run full QIIME workflow generating output files for method comparison
    
    """
    
    # Prepare some variables for the later steps
    fna_fp0 = fna_fps.split(',')[0]
    input_dir, input_filename = split(fna_fp0)
    input_basename, input_ext = splitext(input_filename)
    mapping_categories = parse_mapping_file(open(mapping_fp,'U'))[1]
    if categories:
        # split categories skipping any empty strings (to handle
        # e.g. trailing commas)
        categories = [c for c in categories.split(',') if c]
        for c in categories:
            if c not in mapping_categories:
                raise ValueError, ("Category '%s' is not a column header "
                 "in your mapping file. "
                 "Categories are case and white space sensitive. Valid "
                 "choices are: (%s)" % (c,', '.join(mapping_categories)))
    else:
        categories= []
    if params == None:
        params = parse_qiime_parameters([])
    create_dir(output_dir)
    index_fp = '%s/index.html' % output_dir
    index_links = []
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    log_fp = generate_log_fp(output_dir)
    index_links.append(('Master run log',log_fp,'Log files'))
    logger = WorkflowLogger(log_fp,
                            params=params,
                            qiime_config=qiime_config)
    
    ## Split libraries
    # Prep the split_libraries command
    if suppress_split_libraries:
        if len(fna_fps.split(',')) > 1:
            raise ValueError, \
             "Only one fasta file can be passed when suppress_split_libraries=True"
        split_libraries_seqs_fp = fna_fp0
    else:
        split_libraries_output_dir = \
         '%s/sl_out/' % output_dir
        split_libraries_seqs_fp = \
         '%s/seqs.fna' % split_libraries_output_dir
        split_libraries_hist_fp = \
         '%s/histograms.txt' % split_libraries_output_dir
        split_libraries_log_fp = \
         '%s/split_library_log.txt' % split_libraries_output_dir
        index_links.append(('Demultiplexed sequences',
                            split_libraries_seqs_fp,
                            'Split libraries results'))
        index_links.append(('Split libraries log',
                            split_libraries_log_fp,
                            'Split libraries results'))
        index_links.append(('Sequence length histograms',
                            split_libraries_hist_fp,
                            'Split libraries results'))
        try:
            params_str = get_params_str(params['split_libraries'])
        except KeyError:
            params_str = ''
        # Build the split libraries command
        if qual_fps:
            qual_str = '-q %s' % qual_fps
        else:
            qual_str = ''
        split_libraries_cmd = 'split_libraries.py -f %s %s -m %s -o %s %s' %\
         (fna_fps, qual_str, mapping_fp, split_libraries_output_dir, params_str)
    
        commands.append([('Split libraries', split_libraries_cmd)])
    
        # Call the command handler on the list of commands
        command_handler(commands, 
                        status_update_callback, 
                        logger, 
                        close_logger_on_success=False)
        # Reset the commands list
        commands = []
    
    ## OTU picking through OTU table workflow
    data_analysis_output_dir = '%s/otus/' % output_dir
    de_novo_tree_fp, otu_table_fp = \
     run_qiime_data_preparation(input_fp=split_libraries_seqs_fp, 
                                output_dir=data_analysis_output_dir, 
                                command_handler=command_handler,
                                params=params,
                                qiime_config=qiime_config,
                                parallel=parallel,
                                logger=logger,
                                status_update_callback=status_update_callback)
    index_links.append(('Phylogenetic tree',de_novo_tree_fp,'OTU workflow results'))
    index_links.append(('OTU table',otu_table_fp,'OTU workflow results'))
    
    # If a reference tree was passed, use it for downstream analysis. Otherwise
    # use the de novo tree.
    tree_fp = reference_tree_fp or de_novo_tree_fp
    
    if sampling_depth == None:
        min_count, max_count, median_count, mean_count, counts_per_sample =\
         compute_seqs_per_library_stats(open(otu_table_fp))
        if even_sampling_keeps_all_samples:
            sampling_depth = min_count
        else:
            sampling_depth = \
             guess_even_sampling_depth(counts_per_sample.values())
    
    ## Beta diversity through 3D plots workflow
    bdiv_full_output_dir = '%s/bdiv/' % output_dir
    full_dm_fps = run_beta_diversity_through_plots(
     otu_table_fp=otu_table_fp, 
     mapping_fp=mapping_fp,
     output_dir=bdiv_full_output_dir,
     command_handler=command_handler,
     params=params,
     qiime_config=qiime_config,
     sampling_depth=None,
     tree_fp=tree_fp,
     parallel=parallel,
     logger=logger,
     status_update_callback=status_update_callback)
                            
    # cluster quality stub
    for bdiv_metric, dm_fp in full_dm_fps:
        for category in categories:
            cluster_quality_fp = '%s/%s_%s_cluster_quality.txt' %\
             (bdiv_full_output_dir,bdiv_metric,category)
            try:
                params_str = get_params_str(params['cluster_quality'])
            except KeyError:
                params_str = ''
            # Build the cluster quality command
            cluster_quality_cmd = \
             'cluster_quality.py -i %s -c %s -o %s -m %s %s' %\
              (dm_fp, category, cluster_quality_fp, mapping_fp, params_str)
    
            commands.append([
             ('Cluster quality (%s; %s)' % (bdiv_metric, category),
              cluster_quality_cmd)])
              
            index_links.append(('Cluster quality results (%s, %s)' % (bdiv_metric,category),
                                cluster_quality_fp,
                                'Beta diversity results'))
        # Create links for the bdiv results
        index_links.append(('3D plot (%s, continuous coloring)' % bdiv_metric,
                            '%s/%s_3d_continuous/%s_pc_3D_PCoA_plots.html' % \
                             (bdiv_full_output_dir,bdiv_metric,bdiv_metric),
                            'Beta diversity results'))
        index_links.append(('3D plot (%s, discrete coloring)' % bdiv_metric,
                            '%s/%s_3d_discrete/%s_pc_3D_PCoA_plots.html' % \
                             (bdiv_full_output_dir,bdiv_metric,bdiv_metric),
                            'Beta diversity results'))
        index_links.append(('2D plot (%s, continuous coloring)' % bdiv_metric,
                            '%s/%s_2d_continuous/%s_pc_2D_PCoA_plots.html' % \
                             (bdiv_full_output_dir,bdiv_metric,bdiv_metric),
                            'Beta diversity results'))
        index_links.append(('2D plot (%s, discrete coloring)' % bdiv_metric,
                            '%s/%s_2d_discrete/%s_pc_2D_PCoA_plots.html' % \
                             (bdiv_full_output_dir,bdiv_metric,bdiv_metric),
                            'Beta diversity results'))
        index_links.append(('Distance histograms (%s)' % bdiv_metric,
                            '%s/%s_histograms/%s_dm_distance_histograms.html' % \
                             (bdiv_full_output_dir,bdiv_metric,bdiv_metric),
                            'Beta diversity results'))
        index_links.append(('Distance matrix (%s)' % bdiv_metric,
                            '%s/%s_dm.txt' % \
                             (bdiv_full_output_dir,bdiv_metric),
                            'Beta diversity results'))
        index_links.append(('Principal coordinate matrix (%s)' % bdiv_metric,
                            '%s/%s_pc.txt' % \
                             (bdiv_full_output_dir,bdiv_metric),
                            'Beta diversity results'))

    
    if sampling_depth:
        bdiv_even_output_dir = '%s/bdiv_even%d/' % (output_dir,sampling_depth)
        even_dm_fps = run_beta_diversity_through_plots(
         otu_table_fp=otu_table_fp, 
         mapping_fp=mapping_fp,
         output_dir=bdiv_even_output_dir,
         command_handler=command_handler,
         params=params,
         qiime_config=qiime_config,
         sampling_depth=sampling_depth,
         tree_fp=tree_fp,
         parallel=parallel,
         logger=logger,
         status_update_callback=status_update_callback)
        for bdiv_metric, dm_fp in even_dm_fps:
            for category in categories:
                cluster_quality_fp = '%s/%s_%s_cluster_quality.txt' %\
                 (bdiv_even_output_dir,bdiv_metric,category)
                try:
                    params_str = get_params_str(params['cluster_quality'])
                except KeyError:
                    params_str = ''
                # Build the cluster quality command
                cluster_quality_cmd = \
                 'cluster_quality.py -i %s -c %s -o %s -m %s %s' %\
                  (dm_fp, category, cluster_quality_fp, mapping_fp, params_str)
    
                commands.append([
                 ('Cluster quality (%s; %s)' % (bdiv_metric, category),
                  cluster_quality_cmd)])
                index_links.append(('Cluster quality results (%s, %s)' % (bdiv_metric,category),
                    cluster_quality_fp,
                    'Beta diversity results (even sampling: %d)' % sampling_depth))
                    # Create links for the bdiv results
            index_links.append(('3D plot (%s, continuous coloring)' % bdiv_metric,
                                '%s/%s_3d_continuous/%s_pc_3D_PCoA_plots.html' % \
                                 (bdiv_even_output_dir,bdiv_metric,bdiv_metric),
                                'Beta diversity results (even sampling: %d)' % sampling_depth))
            index_links.append(('3D plot (%s, discrete coloring)' % bdiv_metric,
                                '%s/%s_3d_discrete/%s_pc_3D_PCoA_plots.html' % \
                                 (bdiv_even_output_dir,bdiv_metric,bdiv_metric),
                                'Beta diversity results (even sampling: %d)' % sampling_depth))
            index_links.append(('2D plot (%s, continuous coloring)' % bdiv_metric,
                                '%s/%s_2d_continuous/%s_pc_2D_PCoA_plots.html' % \
                                 (bdiv_even_output_dir,bdiv_metric,bdiv_metric),
                                'Beta diversity results (even sampling: %d)' % sampling_depth))
            index_links.append(('2D plot (%s, discrete coloring)' % bdiv_metric,
                                '%s/%s_2d_discrete/%s_pc_2D_PCoA_plots.html' % \
                                 (bdiv_even_output_dir,bdiv_metric,bdiv_metric),
                                'Beta diversity results (even sampling: %d)' % sampling_depth))
            index_links.append(('Distance histograms (%s)' % bdiv_metric,
                                '%s/%s_histograms/%s_dm_distance_histograms.html' % \
                                 (bdiv_even_output_dir,bdiv_metric,bdiv_metric),
                                'Beta diversity results (even sampling: %d)' % sampling_depth))
            index_links.append(('Distance matrix (%s)' % bdiv_metric,
                                '%s/%s_dm.txt' % \
                                 (bdiv_even_output_dir,bdiv_metric),
                                'Beta diversity results (even sampling: %d)' % sampling_depth))
            index_links.append(('Principal coordinate matrix (%s)' % bdiv_metric,
                                '%s/%s_pc.txt' % \
                                 (bdiv_even_output_dir,bdiv_metric),
                                'Beta diversity results (even sampling: %d)' % sampling_depth))
        
    ## Alpha rarefaction workflow
    arare_full_output_dir = '%s/arare/' % output_dir
    run_qiime_alpha_rarefaction(
     otu_table_fp=otu_table_fp,
     mapping_fp=mapping_fp,
     output_dir=arare_full_output_dir,
     command_handler=command_handler,
     params=params,
     qiime_config=qiime_config,
     tree_fp=tree_fp,
     num_steps=arare_num_steps,
     parallel=parallel,
     logger=logger,
     min_seqs_per_sample=arare_min_seqs_per_sample,
     status_update_callback=status_update_callback)
    
    index_links.append(('Alpha rarefaction plots',
                        '%s/alpha_rarefaction_plots/rarefaction_plots.html'\
                          % arare_full_output_dir,
                        "Alpha rarefaction results"))
    
    
    # OTU category significance and supervised learning
    for category in categories:
        category_signifance_fp = \
         '%s/category_significance_%s.txt' % (output_dir, category)
        try:
            params_str = get_params_str(params['otu_category_significance'])
        except KeyError:
            params_str = ''
        # Build the OTU cateogry significance command
        category_significance_cmd = \
         'otu_category_significance.py -i %s -m %s -c %s -o %s %s' %\
         (otu_table_fp, mapping_fp, category, 
          category_signifance_fp, params_str)
        commands.append([('OTU category significance (%s)' % category, 
                          category_significance_cmd)])
                          
        supervised_learning_dir = \
         '%s/supervised_learning_%s' % (output_dir, category)
        try:
            params_str = get_params_str(params['supervised_learning'])
        except KeyError:
            params_str = ''
        # Build the supervised_learning command
        supervised_learning_cmd = \
         'supervised_learning.py -i %s -m %s -c %s -o %s %s' %\
         (otu_table_fp, mapping_fp, category, 
          supervised_learning_dir, params_str)
        commands.append([('Supervised learning (%s)' % category, 
                          supervised_learning_cmd)])
                          
        index_links.append(('Category significance (%s)' % category,
                    category_signifance_fp,
                    "Category results"))
        index_links.append(('Supervised learning (%s)' % category,
                    supervised_learning_dir,
                    "Supervised learning results"))
    
    command_handler(commands, status_update_callback, logger)
    generate_index_page(index_links,index_fp)

def run_summarize_taxa_through_plots(otu_table_fp, mapping_fp,
    output_dir, mapping_cat, command_handler, params, qiime_config,
    logger=None, status_update_callback=print_to_stdout):
    """ Run the data preparation for summarizing taxonomies and generating plots
    
        The steps performed by this function are:
          1) Summarize OTU by Categore
          2) Summarize Taxonomy
          3) Plot Taxonomy Summary
          
    """
    # Prepare some variables for the later steps
    otu_table_dir, otu_table_filename = split(otu_table_fp)
    otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    create_dir(output_dir)
    
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    
    if logger == None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False
    
    # Prep the summarize otu by category command
    try:
        otu_table_f = open(otu_table_fp,'U')
    except IOError,e:
        logger.write('OTU table filepath cannot be opened. Does it exist?\n' +
                     ' %s\n' % otu_table_fp +
                     'Original Error:\n%s\n' % str(e))
        logger.close()
        raise IOError,e
    
    # if mapping category not passed via command-line, 
    # check if it is passed in params file
    if not mapping_cat:
        try:
            mapping_cat=params['summarize_otu_by_cat']['mapping_category']
        except:
            mapping_cat=None
        
    try:
        params_str = get_params_str(params['summarize_otu_by_cat'])
        # Need to remove the mapping category option, since it is defined above.
        # Using this method since we don't want to change the params dict
        split_params=params_str.split('--')
        updated_params_str=[]
        for i in split_params:
            if not i.startswith('mapping_category'):
                updated_params_str.append(i)
        params_str='--'.join(updated_params_str)
    except:
        params_str = ''
    
    if mapping_cat:
        output_fp=join(output_dir,'%s_otu_table.txt' % (mapping_cat))
        # Build the summarize otu by category command
        summarize_otu_by_cat_cmd = \
         "%s %s/summarize_otu_by_cat.py -i %s -c %s -o %s -m %s %s" %\
         (python_exe_fp, script_dir, mapping_fp, otu_table_fp, output_fp, \
          mapping_cat, params_str)
        
        commands.append(\
         [('Summarize OTU table by Category',summarize_otu_by_cat_cmd)])
         
        otu_table_fp=output_fp
     
    # Prep the summarize taxonomy command
    try:
        params_str = get_params_str(params['summarize_taxa'])
    except:
        params_str = ''
    
    try:
        sum_taxa_levels=params['summarize_taxa']['level']
    except:
        sum_taxa_levels=None
        
    # Build the summarize taxonomy command
    summarize_taxa_cmd = '%s %s/summarize_taxa.py -i %s -o %s %s' %\
     (python_exe_fp, script_dir, otu_table_fp, \
      output_dir, params_str)
    
    commands.append([('Summarize Taxonomy',summarize_taxa_cmd)])

    sum_taxa_fps=[]
    
    if sum_taxa_levels:
        basename=join(output_dir,splitext(split(otu_table_fp)[-1])[0])
        for i in sum_taxa_levels.split(','):
            sum_taxa_fps.append(basename+'_L%s.txt' % (str(i)))
    else:
        basename=join(output_dir,splitext(split(otu_table_fp)[-1])[0])
        # this is the default levels from summarize_taxa, but cannot import
        # script to get these values
        for i in [2,3,4,5,6]:
            sum_taxa_fps.append(basename+'_L%s.txt' % (str(i)))

    # Prep the plot taxa summary plot command(s)
    taxa_summary_plots_dir = '%s/taxa_summary_plots/' % output_dir
    try:
        makedirs(taxa_summary_plots_dir)
    except OSError:
        pass
        
    try:
        params_str = get_params_str(params['plot_taxa_summary'])
    except:
        params_str = ''
    # Build the plot taxa summary plot command(s)

    plot_taxa_summary_cmd =\
         '%s %s/plot_taxa_summary.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, ','.join(sum_taxa_fps),
          taxa_summary_plots_dir, params_str)
    
    commands.append(\
         [('Plot Taxonomy Summary',plot_taxa_summary_cmd)])
    
    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)


def run_ampliconnoise(mapping_fp,
    output_dir, command_handler, params, qiime_config,
    logger=None, status_update_callback=print_to_stdout,
    chimera_alpha=-3.8228,chimera_beta=0.6200, sff_txt_fp=None, numnodes=2,
    suppress_perseus=True, output_filepath=None, platform='flx',
    seqnoise_resolution=None, truncate_len=None):
    """ Run the ampliconnoise pipeline
    
        The steps performed by this function are:
1. Split input sff.txt file into one file per sample

2. Run scripts required for PyroNoise

3. Run scripts required for SeqNoise

4. Run scripts requred for Perseus (chimera removal)

5. Merge output files into one file similar to the output of split_libraries.py

    output_filepath should be absolute
    seqnoise_resolution should be string
    environment variable PYRO_LOOKUP_FILE must be set correctly. Thus be
    careful passing command handlers that don't spawn child processes, as they
    may not inherit the correct environment variable setting
    """
    map_data,headers,comments = parse_mapping_file(open(mapping_fp,'U'))
    create_dir(output_dir)

    if seqnoise_resolution == None:
        if platform=='flx': seqnoise_resolution = '30.0'
        elif platform=='titanium': seqnoise_resolution = '25.0'
        else: raise RuntimeError('seqnoise_resolution not set, and no'+\
            ' default for platform '+platform)

    if truncate_len == None:
        if platform=='flx': truncate_len = '220'
        elif platform=='titanium': truncate_len = '400'
        else: raise RuntimeError('truncate_len not set, and no'+\
            ' default for platform '+platform)

    sample_names = [] # these are filenames minus extension, and are sample IDs
    primer_seqs = [] # same order as sample_names
    bc_seqs = [] # same order as sample_names
    for i in range(len(map_data)):
        sample_names.append(map_data[i][headers.index('SampleID')])
        bc_seqs.append(map_data[i][headers.index('BarcodeSequence')])
        # don't know why don't just take off the primer now. 
        # but that's done later
        # primer += (map_data[i][headers.index('LinkerPrimerSequence')])
        # for char, bases in IUPAC_DNA_ambiguities.items():
        #     primer = primer.replace(char,'['+''.join(bases)+']')

        primer = (map_data[i][headers.index('LinkerPrimerSequence')])
        for char, bases in IUPAC_DNA_ambiguities.items():
            primer = primer.replace(char,'['+''.join(bases)+']')
        primer_seqs.append(primer)

    if len(set(primer_seqs)) != 1:
        raise RuntimeError(
            'Error: only one primer per mapping file supported.')
    one_primer = primer_seqs[0]

    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()

    if logger == None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False


    # execute commands in output_dir
    called_dir = os.getcwd()
    os.chdir(output_dir)
    fh = open(os.path.join(output_dir,'map.csv'),'w')
    for i in range(len(sample_names)):
        fh.write(sample_names[i]+','+bc_seqs[i]+'\n')
    fh.close()

    # these are the fasta results, e.g. PC.636_Good.fa
    # later we merge them and copy to output file
    post_pyro_tail = '_'+truncate_len
    if suppress_perseus == True:
        fasta_result_names = [sample_name + post_pyro_tail+'_seqnoise_cd.fa'
          for sample_name in sample_names]
    else:
        fasta_result_names = [sample_name + '_Good.fa' \
          for sample_name in sample_names]

    cmd = 'cd '+output_dir # see also os.chdir above
    commands.append([('change to output dir', cmd)])
    
    cmd = 'echo $PYRO_LOOKUP_FILE > pyro_lookup_filepath.txt'
    commands.append([('confirm pyro lookup filepath environment variable',
        cmd)])


    cmd = 'SplitKeys.pl '+one_primer+' map.csv < '+\
        os.path.join(called_dir,sff_txt_fp)+\
        ' > splitkeys_log.txt 2> unassigned.fna'
    commands.append([('split sff.txt via barcodes (keys)', cmd)])

    for i, sample_name in enumerate(sample_names):

        # Build the summarize taxonomy command
        if platform == 'flx':
            cmd = 'Clean360.pl '+one_primer+' '+sample_name+' < '+\
                sample_name+'.raw'
            commands.append([('clean flows '+sample_name, cmd)])

            # these run through the whole sff file once per sample, I think
            # cmd = "FlowsFA.pl " + primer_seqs[i] + ' '+sample_name +' < '+\
            #     os.path.join(called_dir,sff_txt_fp)
            # commands.append([('extract flows '+sample_name, cmd)])
        elif platform == 'titanium':
            cmd = 'CleanMinMax.pl '+one_primer+' '+sample_name+' < '+\
                sample_name+'.raw'
            commands.append([('clean flows '+sample_name, cmd)])

            # cmd = "FlowsMinMax.pl " + primer_seqs[i] + ' '+sample_name +' < '+\
            #     os.path.join(called_dir,sff_txt_fp)
            # commands.append([('extract flows '+sample_name, cmd)])
        else:
            raise RuntimeError("platform " + platform + " not supported")

        cmd = "mpirun -np "+str(numnodes)+" PyroDist -in "+\
          sample_name+".dat -out "+sample_name+ " > "+sample_name+".pdout"
        commands.append([('pyrodist '+sample_name, cmd)])

        cmd = "FCluster -in "+sample_name+".fdist -out "+sample_name+\
          " > "+sample_name+".fcout"
        commands.append([('fcluster pyrodist '+sample_name, cmd)])

# e.g.:
# mpirun -np 2 PyroNoise -din PC.354.dat -out PC.354_pyronoise -lin
# PC.354.list -s 60.0 -c 0.01 > PC.354_pyronoise.pnout
        cmd = "mpirun -np "+str(numnodes)+" PyroNoise -din "+\
            sample_name+".dat -out "+\
            sample_name+"_pyronoise "+"-lin "+\
            sample_name+".list -s 60.0 -c 0.01 > "+\
            sample_name+"_pyronoise.pnout"
        commands.append([('pyronoise '+sample_name, cmd)])

        cmd = 'Parse.pl '+bc_seqs[i]+' '+truncate_len+' < '+\
            sample_name+'_pyronoise_cd.fa'+' > '+ sample_name+'_'+\
            truncate_len+'.fa'
        commands.append([('truncate '+sample_name, cmd)])

        # now start with post_pyro_tail
        cmd = "mpirun -np "+str(numnodes)+" SeqDist -in "+\
            sample_name+post_pyro_tail+\
            ".fa > "+sample_name+post_pyro_tail+".seqdist"
        commands.append([('seqdist '+sample_name, cmd)])

        cmd = "FCluster -in "+sample_name+post_pyro_tail+".seqdist -out "+\
            sample_name+post_pyro_tail+"fcl > "+\
            sample_name+post_pyro_tail+".fcout"
        commands.append([('fcluster seqdist '+sample_name, cmd)])

# e.g.:
# mpirun -np 2 SeqNoise -in PC.354_pyronoise_cd.fa -din
# PC.354_pyronoise_cd.seqdist -out PC.354_pyronoise_cd_seqnoise -lin
# PC.354_pyronoise_cdfcl.list -min PC.354_pyronoise.mapping -s 30.0 -c 0.08 >
# PC.354_pyronoise_cd.snout

        cmd = "mpirun -np "+str(numnodes)+" SeqNoise -in "+\
            sample_name+post_pyro_tail+\
            ".fa -din "+sample_name+post_pyro_tail+".seqdist -out "+\
            sample_name+post_pyro_tail+\
            "_seqnoise -lin "+sample_name+post_pyro_tail+'fcl.list -min '+\
            sample_name+'_pyronoise'+\
            '.mapping -s '+seqnoise_resolution+' -c 0.08 > '+\
            sample_name+post_pyro_tail+'.snout'
        commands.append([('seqnoise '+sample_name, cmd)])


        if suppress_perseus == False: 

            cmd = 'Perseus -sin '+sample_name+post_pyro_tail+\
                '_seqnoise_cd.fa > ' +\
                sample_name+'.per'
            commands.append([('Perseus '+sample_name, cmd)])

            cmd = 'Class.pl '+sample_name+'.per '+\
                str(chimera_alpha) + ' '+ str(chimera_beta)+\
                ' > '+sample_name+'.class'
            commands.append([('Class.pl '+sample_name, cmd)])

            cmd = 'FilterGoodClass.pl '+sample_name+post_pyro_tail+\
                '_seqnoise_cd.fa '+\
                sample_name+'.class 0.5 > '+sample_name+'_Chi.fa 2> '+\
                sample_name+'_Good.fa'
            commands.append([('FilterGoodClass '+sample_name, cmd)])

        cmd = '%s %s/unweight_fasta.py -i %s -o %s -l %s' %\
         (python_exe_fp, script_dir, fasta_result_names[i], 
         sample_name+'_unw.fna', sample_name)
        commands.append([('unweight fasta '+sample_name, cmd)])

    cmd = 'cat ' +\
      ' '.join([sample_name+'_unw.fna' for sample_name in sample_names]) +\
      ' > ' + output_filepath # this should be an abs filepath
    commands.append([('cat into one fasta file', cmd)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)

