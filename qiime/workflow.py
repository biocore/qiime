#!/usr/bin/env python
# File created on 30 Dec 2009.
from __future__ import division
from subprocess import Popen, PIPE, STDOUT
from os import makedirs, listdir
from glob import glob
from os.path import split, splitext, join, dirname, abspath
from datetime import datetime
from qiime.parse import parse_mapping_file
from qiime.util import (compute_seqs_per_library_stats, 
                        get_qiime_scripts_dir,
                        create_dir)
from qiime.make_sra_submission import twocol_data_to_dict, read_tabular_data, generate_output_fp
from qiime.sra_spreadsheet_to_map_files import get_study_groups

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

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
        self.write('Logging started at %s\n\n' % start_time)
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
                   logger):
    """Print list of commands to run """
    logger.write("Printing commands only.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('#%s' % e[0])
            print '%s' % e[1]
            logger.write('# %s command\n%s\n\n' % e)
            
def call_commands_serially(commands,
                           status_update_callback,
                           logger):
    """Run list of commands, one after another """
    logger.write("Executing commands.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('%s\n%s' % e)
            logger.write('# %s command \n%s\n\n' % e)
            proc = Popen(e[1],shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
            return_value = proc.wait()
            if return_value != 0:
                msg = "\n\n*** ERROR RAISED DURING STEP: %s\n" % e[0] +\
                 "Command run was:\n %s\n" % e[1] +\
                 "Command returned exit status: %d\n" % return_value +\
                 "Stdout/stderr:\n%s\n" % proc.stdout.read()
                logger.write(msg)
                logger.close()
                raise WorkflowError, msg
    logger.close()

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

## End utilities used by the workflow functions

## Begin task-specific workflow functions
def run_qiime_data_preparation(input_fp, output_dir, command_handler,
    params, qiime_config, sff_input_fp=None, mapping_fp=None,
    parallel=False, status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime 
    
        The steps performed by this function are:
          0) Optionally denoise the sequences (if sff_input_fp=True);
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
    logger = WorkflowLogger(generate_log_fp(output_dir),
                            params=params,
                            qiime_config=qiime_config)
    
    # Prep the denoising command
    if sff_input_fp != None:
        denoise = True
        assert mapping_fp != None,\
         "Mapping file must be provided for denoising."+\
         " (Need to extract the primer sequence.)"
        denoise_output_dir = '%s/denoised_seqs/' % output_dir
        denoised_seqs_fp = '%s/denoised_seqs.fasta' % denoise_output_dir
        denoised_mapping_fp = '%s/denoiser_mapping.txt' % denoise_output_dir
        
        if parallel:
            parallel_str = '-n %s' % qiime_config['jobs_to_start']
        else:
            parallel_str = ''
            
        try:
            params_str = get_params_str(params['denoise'])
        except KeyError:
            params_str = ''
        
        # build the denoiser command
        denoise_cmd = '%s %s/denoise.py -i %s -f %s --method fast -m %s -o %s %s %s' %\
         (python_exe_fp, script_dir, sff_input_fp, input_fp, mapping_fp,
          denoise_output_dir, parallel_str, params_str)
        commands.append([('Denoise', denoise_cmd)])
        
        # some values that get passed to subsequent steps change when 
        # denoising -- set those here
        original_input_fp = input_fp
        input_fp = denoised_seqs_fp
        input_basename, input_ext = splitext(split(denoised_seqs_fp)[1])
    else:
        denoise = False
    
    # Prep the OTU picking command
    otu_picking_method = params['pick_otus']['otu_picking_method']
    pick_otu_dir = '%s/%s_picked_otus' % (output_dir, otu_picking_method)
    otu_fp = '%s/%s_otus.txt' % (pick_otu_dir,input_basename)
    if parallel and otu_picking_method == 'blast':
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
            
        # Build the OTU picking command
        pick_otus_cmd = '%s %s/parallel_pick_otus_blast.py -i %s -o %s -T %s' %\
         (python_exe_fp, script_dir, input_fp, pick_otu_dir, params_str)
    else:
        if denoise:
            # we want to make sure the user is using the right set of commands
            # For now we force to use uclust --user_sort --optimal
            # in the future we might want to do this more clever
            # and force the user to have a good parameter set in the config file
            if 'optimal_uclust' not in params['pick_otus']:
                logger.write("Warning: Setting option pick_otus:optimal_uclust to True "
                             + "for compatibility with denoising\n")
            params['pick_otus']['optimal_uclust']=None

            if 'user_sort' not in params['pick_otus']:
                logger.write("Warning: Setting option pick_otus:user_sort to True "
                                 + "for compatibility with denoising\n")
            params['pick_otus']['user_sort']=None

            if 'presort_by_abundance_uclust' in params['pick_otus']:
                logger.write("Warning: Disabling option pick_otus:presort_by_abundance_uclust "
                              +"with uclust OTU picker for compatibility with denoising")
                del params['pick_otus']['presort_by_abundance_uclust']
        try:
            params_str = get_params_str(params['pick_otus'])
        except KeyError:
            params_str = ''
        # Build the OTU picking command
        pick_otus_cmd = '%s %s/pick_otus.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, input_fp, pick_otu_dir, params_str)

    commands.append([('Pick OTUs', pick_otus_cmd)])
    
    # Prep the merge_denoiser_output.py command, if denoising
    if denoise:
        pick_otu_dir = '%s/denoised_otus/' % pick_otu_dir
        
        try:
            params_str = get_params_str(params['merge_denoiser_output'])
        except KeyError:
            params_str = ''
        merge_denoiser_output_cmd = \
         '%s %s/merge_denoiser_output.py -m %s -p %s -f %s -d %s -o %s %s' %\
         (python_exe_fp, script_dir, denoised_mapping_fp, otu_fp, 
          original_input_fp, denoised_seqs_fp, pick_otu_dir, params_str)
          
        input_fp = '%s/denoised_all.fasta' % pick_otu_dir
        otu_fp = '%s/denoised_otu_map.txt' % pick_otu_dir
        commands.append([('Merge denoiser output', merge_denoiser_output_cmd)])
    
    # Prep the representative set picking command
    rep_set_dir = '%s/rep_set/' % pick_otu_dir
    try:
        makedirs(rep_set_dir)
    except OSError:
        pass
    rep_set_fp = '%s/%s_rep_set.fasta' % (rep_set_dir,input_basename)
    rep_set_log_fp = '%s/%s_rep_set.log' % (rep_set_dir,input_basename)
    
    if denoise:
        #force rep_set picking methd to be 'first' if not already set
        #Required for picking output from merge_denoiser_output
        if ('rep_set_picking_method' in params['pick_rep_set']
            and not params['pick_rep_set']['rep_set_picking_method'] == 'first'):
            logger.write("Warning: Setting pick_rep_set:rep_set_picking_method to 'first' "+
                     "for compatibility with denoising.\n")
            params['pick_rep_set']['rep_set_picking_method'] = 'first'
        
    try:
        params_str = get_params_str(params['pick_rep_set'])
    except KeyError:
        params_str = ''
    # Build the representative set picking command
    pick_rep_set_cmd = '%s %s/pick_rep_set.py -i %s -f %s -l %s -o %s %s' %\
     (python_exe_fp, script_dir, otu_fp, input_fp, rep_set_log_fp,\
      rep_set_fp, params_str)
    commands.append([('Pick representative set', pick_rep_set_cmd)])
    
    # Prep the pynast alignment command
    pynast_dir = '%s/%s_aligned_seqs' % \
     (rep_set_dir,params['align_seqs']['alignment_method'])
    aln_fp = '%s/%s_rep_set_aligned.fasta' % (pynast_dir,input_basename)
    alignment_method = params['align_seqs']['alignment_method']
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

    
    # Prep the taxonomy assignment command
    assignment_method = params['assign_taxonomy']['assignment_method']
    assign_taxonomy_dir = '%s/%s_assigned_taxonomy' %\
     (rep_set_dir,assignment_method)
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
    
    # Append commands which can be run simulataneously in parallel
    commands.append([('Align sequences', align_seqs_cmd),\
                     ('Assign taxonomy',assign_taxonomy_cmd)])
    
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
     (python_exe_fp, script_dir, filtered_aln_fp, tree_fp, log_fp,\
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
     (python_exe_fp, script_dir, otu_fp, taxonomy_fp, otu_table_fp, params_str)
    
    # Append commands which can be run simulataneously in parallel
    commands.append([('Build phylogenetic tree', make_phylogeny_cmd),\
                     ('Make OTU table', make_otu_table_cmd)])
    
    # Call the command handler on the list of commands
    command_handler(commands,status_update_callback,logger=logger)
    
def run_beta_diversity_through_3d_plot(otu_table_fp, mapping_fp,\
    output_dir, command_handler, params, qiime_config, tree_fp=None,\
    parallel=False, status_update_callback=print_to_stdout):
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
    logger = WorkflowLogger(generate_log_fp(output_dir),
                            params=params,
                            qiime_config=qiime_config)
    
    mapping_file_header = parse_mapping_file(open(mapping_fp,'U'))[1]
    mapping_fields = ','.join(mapping_file_header)
    
    beta_diversity_metrics = params['beta_diversity']['metrics'].split(',')
    
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
    
    # Prep the 3d prefs file generator command
    prefs_fp = '%s/prefs.txt' % output_dir
    try:
        params_str = get_params_str(params['make_prefs_file'])
    except KeyError:
        params_str = ''
    # Build the 3d prefs file generator command
    prefs_cmd = \
     '%s %s/make_prefs_file.py -m %s -o %s %s' %\
     (python_exe_fp, script_dir, mapping_fp, prefs_fp, params_str)
    commands.append([('Build prefs file', prefs_cmd)])
        
    for beta_diversity_metric in beta_diversity_metrics:
        
        beta_div_fp = '%s/%s_%s' % \
         (output_dir, beta_diversity_metric, otu_table_filename)
        
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
    
    # Call the command handler on the list of commands
    command_handler(commands, status_update_callback, logger)


def run_qiime_alpha_rarefaction(otu_table_fp, mapping_fp,\
    output_dir, command_handler, params, qiime_config, tree_fp=None,\
    num_steps=10, parallel=False, min_seqs_per_sample=10,\
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
    logger = WorkflowLogger(generate_log_fp(output_dir),
                            params=params,
                            qiime_config=qiime_config)
    
    alpha_diversity_metrics = params['alpha_diversity']['metrics'].split(',')
    
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
    step = int((median_count - min_seqs_per_sample) / num_steps)
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
        try:
            # Want to find a cleaner strategy for this: the rarefaction 
            # parallel script doesn't support the jobs_to_start option -
            # one job is started per rarefied otu table to be created -
            # so need to remove this option. This works for now though.
            d = params['parallel'].copy()
            del d['jobs_to_start']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass        
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
        try:
            # Want to find a cleaner strategy for this: the alpha diversity 
            # parallel script doesn't support the jobs_to_start option -
            # one job is started per rarefied otu table to be created -
            # so need to remove this option. This works for now though.
            d = params['parallel'].copy()
            del d['jobs_to_start']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass   
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
    command_handler(commands,status_update_callback,logger)

    
def run_jackknifed_upgma_clustering(otu_table_fp,tree_fp,seqs_per_sample,\
    output_dir, command_handler, params, qiime_config,\
    parallel=False,status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime 
    
        The steps performed by this function are:
          1) Compute beta diversity distance matrix from otu table (and
           tree, if applicable)
          2) Build rarefied OTU tables;
          3) Build UPGMA tree from full distance matrix;
          4) Compute distance matrics for rarefied OTU tables;
          5) Build UPGMA trees from rarefied OTU table distance matrices;
          6) Compare rarefied OTU table distance matrix UPGMA trees 
           to tree full UPGMA tree and write support file and newick tree
           with support values as node labels.
    """
    # Prepare some variables for the later steps
    otu_table_dir, otu_table_filename = split(otu_table_fp)
    otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    create_dir(output_dir)
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    logger = WorkflowLogger(generate_log_fp(output_dir),
                            params=params,
                            qiime_config=qiime_config)
    
    beta_diversity_metrics = params['beta_diversity']['metrics'].split(',')
    
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
    #     try:
    #         # Want to find a cleaner strategy for this: the rarefaction 
    #         # parallel script doesn't support the jobs_to_start option -
    #         # one job is started per rarefied otu table to be created -
    #         # so need to remove this option. This works for now though.
    #         d = params['parallel'].copy()
    #         del d['jobs_to_start']
    #         params_str += ' %s' % get_params_str(d)
    #     except KeyError:
    #         pass        
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
        master_tree_fp = '%s/%s_upgma.tre' % (metric_output_dir,otu_table_basename)
        try:
            params_str = get_params_str(params['upgma_cluster'])
        except KeyError:
            params_str = ''
        # Build the hierarchical clustering command (for full distance matrix)
        hierarchical_cluster_cmd = '%s %s/upgma_cluster.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, distance_matrix_fp, master_tree_fp, params_str)
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
            try:
                # Want to find a cleaner strategy for this: the beta diversity 
                # parallel script doesn't support the jobs_to_start option -
                # one job is started per rarefied otu table to be created -
                # so need to remove this option. This works for now though.
                d = params['parallel'].copy()
                del d['jobs_to_start']
                params_str += ' %s' % get_params_str(d)
            except KeyError:
                pass        
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
        tree_compare_cmd = '%s %s/tree_compare.py -s %s -m %s -o %s %s' %\
         (python_exe_fp, script_dir, upgma_dir, master_tree_fp, \
          tree_compare_dir, params_str)
        commands.append(\
         [('Tree compare (%s)' % beta_diversity_metric,\
           tree_compare_cmd)])
           
    # Call the command handler on the list of commands
    command_handler(commands,status_update_callback,logger)


## Begin SRA submission workflow and related functions

def get_submission_info(submission_fp):
    f = open(submission_fp, 'U')
    return twocol_data_to_dict(read_tabular_data(f)[1])

def get_run_info(experiment_fp):
    infile = open(experiment_fp, 'U')
    _, study_groups = get_study_groups(infile)
    return study_groups.keys()            

def get_sff_filenames(sff_dir, run_prefix):
    return filter(
     lambda x: x.startswith(run_prefix) and x.endswith('.sff'),
     listdir(sff_dir))



def run_process_sra_submission(
    input_experiment_fp,
    input_submission_fp,
    sff_dir,
    refseqs_fp,
    output_dir,
    params,
    qiime_config,
    command_handler,
    status_update_callback=print_to_stdout,
    remove_unassigned=[],
    experiment_link_fp=None,
    experiment_attribute_fp=None):
    """Run the SRA second-stage submission process.

    If human screening should be bypassed, pass refseqs_fp=None.

    The steps performed by this function are:
        Get fasta and qual from sff files
        Produce valid mapping file for library demultiplexing
        Demultiplex libraries
        Optional human screen: Pick otus with uclust_ref against 
           reference database, discarding sequences that don't 
           hit the reference database. The resulting otu map is then 
           the sequence identifiers which pass the human screen.
        Make per-library files of good ids to pass to sfffile
        Use sfffile to make per-library sff files
        Use sfffile to quality-trim the barcodes, primers and linkers
        Move files around and make archive
        Finally, make the XML files for a second-stage submission

    The arguments input_experiment_fp, input_submission_fp, sff_dir,
     experiment_link_fp, and experiment_attribute_fp have corresponding
     arguments in make_sra_submission.py.

    The refseqs_fp argument corresponds to the querydb argument of
     exclude_seqs_by_blast.py.  It is to be the path to a FASTA file of 
     reference sequences. If human screening should be bypassed, 
     pass refseqs_fp=None.

    The remove_unassigned keyword argument is a list of run prefixes
     for which to remove unassigned sequences.
    """
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    submission_dir = dirname(input_experiment_fp)
    create_dir(output_dir)
    # update log_fp to go into the directory that
    # makes the most sense 
    logger = WorkflowLogger(generate_log_fp(output_dir),
                            params=params,
                            qiime_config=qiime_config)

    submission_info = get_submission_info(input_submission_fp)
    if 'file' in submission_info:
        # if a sff tar filename was provided in the submission,
        # grab it and copy the submission file
        submission_tar_fn = submission_info['file']
        second_stage_submission_fp = join(output_dir,split(input_submission_fp)[1])
        commands.append([(
            'Create a copy of submission text file in output directory',
            'cp %s %s' % (input_submission_fp, second_stage_submission_fp))])
    else:
        # if a sff tar filename was not provided in the submission, create
        # a name from the submission_id, and append it to the copy of the
        # submission file
        submission_tar_fn = \
         '%s.tgz' % submission_info['submission_id'].replace(' ','_')
        second_stage_submission_fp = join(output_dir,'submission_second_stage.txt')
        second_stage_submission_f = open(second_stage_submission_fp,'w')
        second_stage_submission_f.write(open(input_submission_fp,'U').read())
        second_stage_submission_f.write('\n%s' %
         '\t'.join(['file',submission_tar_fn,submission_tar_fn,
                    "tgz filename, if submitting sffs"]))
        second_stage_submission_f.close()
    
    submission_tar_fp = join(output_dir, submission_tar_fn)

    # Prelude: Create sff directory for submission data
    submission_sff_dir = join(output_dir, 'per_run_sff')
    create_dir(submission_sff_dir)
    sff_working_dir = join(output_dir, 'sff_files')
    create_dir(sff_working_dir)
    
    # why are we creating copies of the input file?
    input_experiment_copy_fp = generate_output_fp(
        input_experiment_fp, '.txt', output_dir)
    commands.append([(
        'Create a copy of experiment text file in output directory',
        'cp %s %s' % (input_experiment_fp, input_experiment_copy_fp))])
        
    if refseqs_fp:
        refseqs_copy_fp = generate_output_fp(refseqs_fp, '.fasta', output_dir)
        commands.append([(
            'Create a copy of reference set FASTA file in output directory',
            'cp %s %s' % (refseqs_fp, refseqs_copy_fp))])
        perform_human_screen = True
        logger.write('Human screening will be performed against %s.\n' % refseqs_fp)
    else:
        perform_human_screen = False
        logger.write('No reference sequences provided.'
                     ' HUMAN SCREENING WILL NOT BE PERFORMED!\n')
        
        
    if not sff_dir.endswith('/'):
        sff_dir = sff_dir + '/'
    bash_command = '"cp %s*.sff %s"' % (sff_dir, sff_working_dir)
    # why do we need to call bash here?
    commands.append([(
        'Create a copy of sff files in working directory',
        ('bash -c %s' % bash_command))])

    # Step 1
    process_sff_cmd = '%s %s/process_sff.py -i %s' %\
     (python_exe_fp,script_dir, sff_working_dir)

    commands.append([(
        'Process SFF files to create FASTA and QUAL files',
        process_sff_cmd)])

    # Step 2
    sra_spreadsheet_to_map_files_cmd = \
     '%s %s/sra_spreadsheet_to_map_files.py -i %s' %\
     (python_exe_fp,script_dir, input_experiment_copy_fp)
    commands.append([(
        'Create mapping files from the SRA experiment input file',
        sra_spreadsheet_to_map_files_cmd)])

    for study_ref, run_prefix in get_run_info(input_experiment_fp):

        # Step 3: split libraries
        map_fp = join(
            output_dir, '%s_%s.map' % (study_ref, run_prefix))
        sff_filenames = get_sff_filenames(sff_dir, run_prefix)
        sff_basenames = [splitext(x)[0] for x in sff_filenames]
        sff_basepaths = [join(sff_working_dir, x) for x in sff_basenames]
        fna_string = ','.join([b + '.fna' for b in sff_basepaths])
        qual_string = ','.join([b + '.qual' for b in sff_basepaths])
        library_dir = join(output_dir, '%s_demultiplex' % run_prefix)
        
        params_str = get_params_str(params['split_libraries'])
        
        split_libraries_cmd = \
         '%s %s/split_libraries.py -m %s -f %s -q %s -o %s %s' %\
         (python_exe_fp,script_dir,map_fp,fna_string, qual_string,
         library_dir, params_str)
        
        ## WHY IS THIS BEING DONE? SHOULD THIS JUST BE ADDED TO THE 
        ## PARAMETERS FILE? PROBABLY...
        if run_prefix in remove_unassigned:
            split_libraries_args.append('-r')
        
        commands.append([(
            'Demultiplex run %s' % run_prefix, split_libraries_cmd)])
        seqs_fp = join(library_dir, 'seqs.fna')

        if perform_human_screen:
            # pick_otus against reference set for human screen
            params_str = get_params_str(params['pick_otus'])
            pick_otu_params = set(params['pick_otus'].keys())
            if pick_otu_params and\
               pick_otu_params != set(['enable_rev_strand_match','similarity']):
                raise WorkflowError,\
                 ("pick_otus only supports passing of similarity and"
                  " enable_rev_strand_match during SRA submission workflow.")
        
            pick_otus_cmd = \
             '%s %s/pick_otus.py -i %s -o %s -m uclust_ref --suppress_new_clusters -r %s --max_accepts 1 --suppress_presort_by_abundance_uclust --user_sort %s' %\
             (python_exe_fp, script_dir, seqs_fp, library_dir, refseqs_copy_fp, params_str)
            commands.append([('Human screen with uclust_ref OTU picker', pick_otus_cmd)])

            # Step xx - Screen input seqs to filter human sequences
            otus_fp = join(library_dir, 'seqs_otus.txt')        
            screened_seqs_fp = join(library_dir, 'screened_seqs.fasta')
            params_str = get_params_str(params['filter_fasta'])
            filter_fasta_cmd = \
             '%s %s/filter_fasta.py -f %s -o %s -m %s %s' % \
             (python_exe_fp, script_dir, seqs_fp, screened_seqs_fp, otus_fp, params_str)

            commands.append([(
                'Filter input sequences to remove those which didn\'t pass human screen',
                filter_fasta_cmd)])
        else:
            screened_seqs_fp = seqs_fp
        
        # Step 7 - make per library id lists
        per_lib_sff_dir = join(library_dir, 'per_lib_info')
        params_str = get_params_str(params['make_library_id_lists'])
        make_library_id_lists_cmd = \
         '%s %s/make_library_id_lists.py -i %s -o %s %s' % \
         (python_exe_fp, script_dir, screened_seqs_fp, per_lib_sff_dir, params_str)

        commands.append([(
            'Create per-library id lists to use when splitting SFF files',
            make_library_id_lists_cmd)])

        # Step 8 -- make per library sff files
        sff_string = ','.join(
                    [join(sff_working_dir, x) for x in sff_filenames])
        params_str = get_params_str(params['make_per_library_sff'])
        make_per_library_sff_cmd = \
         '%s %s/make_per_library_sff.py -i %s -l %s %s' %\
         (python_exe_fp, script_dir, sff_string, per_lib_sff_dir, params_str)
        commands.append([(
            'Create per-library SFF files', make_per_library_sff_cmd)])

        # Step 9 -- trim sff primers
        params_str = get_params_str(params['trim_sff_primers'])
        trim_sff_primers_cmd = \
         '%s %s/trim_sff_primers.py -m %s -l %s %s' %\
         (python_exe_fp, script_dir, map_fp, per_lib_sff_dir, params_str)

        commands.append([(
            'Trim primer sequences from per-library SFF files',
            trim_sff_primers_cmd)])

        # Step 10 -- organize submission files

        run_sff_output_dir = join(submission_sff_dir, run_prefix)
        create_dir(run_sff_output_dir)

        # for sff_fp in glob(join(per_lib_sff_dir, '*.sff')):
        #     filename = split(fp)[1]
        #     rename(fp,'%s/%s' % (run_sff_submission_dir,filename))
        if not per_lib_sff_dir.endswith('/'):
            per_lib_sff_dir = per_lib_sff_dir + '/'
        if not run_sff_output_dir.endswith('/'):
            run_sff_output_dir = run_sff_output_dir + '/'
        bash_command = (
            '"cp %s*.sff %s"' % (per_lib_sff_dir, run_sff_output_dir))
        commands.append([(
             'Copy per-library SFF files to submission directory',
            'bash -c %s' % bash_command)])

        orig_unassigned_fp = join(run_sff_output_dir, 'Unassigned.sff')
        desired_unassigned_fp = join(
            run_sff_output_dir, '%s_default_%s.sff' % (study_ref, run_prefix))
        commands.append([('Rename Unassigned.sff',
         'mv %s %s' % (orig_unassigned_fp, desired_unassigned_fp))])
    
    commands.append([('Create archive of per-library SFF files',
                     'cd %s ; tar -czf %s %s' % 
                     (abspath(output_dir), split(submission_tar_fp)[1], 
                      split(submission_sff_dir)[1]))])

    # Step 11
    params_str = get_params_str(params['make_sra_submission'])
    make_sra_submission_cmd = \
     '%s %s/make_sra_submission.py -u %s -e %s -s %s -o %s %s' %\
     (python_exe_fp, script_dir, second_stage_submission_fp, input_experiment_copy_fp,
      submission_sff_dir, output_dir, params_str)
    commands.append([('Make SRA submission XML files', 
                      make_sra_submission_cmd)])

    # Call the command handler on the list of commands
    command_handler(commands, status_update_callback, logger)

    
## End task-specific workflow functions
    
