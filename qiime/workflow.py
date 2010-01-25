#!/usr/bin/env python
# File created on 30 Dec 2009.
from __future__ import division
from subprocess import Popen, PIPE, STDOUT
from os import makedirs
from glob import glob
from os.path import split, splitext
from qiime.parse import parse_map
from qiime.util import compute_seqs_per_library_stats

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__status__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

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
            proc = Popen(e[1],shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
            return_value = proc.wait()
            if return_value != 0:
                msg = "\n\n*** ERROR RAISED DURING STEP: %s\n" % e[0] +\
                 "Command run was:\n %s\n" % e[1] +\
                 "Command returned exit status: %d\n" % return_value +\
                 "Stdout/stderr:\n%s\n" % proc.stdout.read()
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
        pick_otus_cmd = '%s %s/parallel/pick_otus_blast.py -i %s -o %s -T %s' %\
         (python_exe_fp, qiime_dir, input_fp, pick_otu_dir, params_str)
    else:
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
        align_seqs_cmd = '%s %s/parallel/align_seqs_pynast.py -i %s -o %s -T %s' %\
         (python_exe_fp, qiime_dir, rep_set_fp, pynast_dir, params_str)
    else:
        try:
            params_str = get_params_str(params['align_seqs'])
        except KeyError:
            params_str = ''
        # Build the pynast alignment command
        align_seqs_cmd = '%s %s/align_seqs.py -i %s -o %s %s' %\
         (python_exe_fp, qiime_dir, rep_set_fp, pynast_dir, params_str)

    
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
         '%s %s/parallel/assign_taxonomy_%s.py -i %s -o %s -T %s' %\
         (python_exe_fp, qiime_dir, assignment_method, rep_set_fp,\
          assign_taxonomy_dir, params_str)
    else:
        try:
            params_str = get_params_str(params['assign_taxonomy'])
        except KeyError:
            params_str = ''
        # Build the taxonomy assignment command
        assign_taxonomy_cmd = '%s %s/assign_taxonomy.py -o %s -i %s %s' %\
         (python_exe_fp, qiime_dir, assign_taxonomy_dir,\
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
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    qiime_home = qiime_config['qiime_home']
    qiime_dir = qiime_config['qiime_dir']
    
    mapping_file_header = parse_map(open(mapping_fp,'U'),return_header=True)[0][0]
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
     (python_exe_fp, qiime_dir, otu_table_fp, output_dir, params_str)
    commands.append(\
     [('Beta Diversity (%s)' % ', '.join(beta_diversity_metrics), beta_div_cmd)])
    
    # Prep the 3d prefs file generator command
    prefs_fp = '%s/3d_prefs.txt' % output_dir
    try:
        params_str = get_params_str(params['make_3d_plot_prefs_file'])
    except KeyError:
        params_str = ''
    # Build the 3d prefs file generator command
    prefs_cmd = \
     '%s %s/../scripts/make_3d_plot_prefs_file.py -b "%s" -p %s %s' %\
     (python_exe_fp, qiime_dir, mapping_fields, prefs_fp, params_str)
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
         (python_exe_fp, qiime_dir, beta_div_fp, pc_fp, params_str)
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
          (python_exe_fp, qiime_dir, prefs_fp, pc_fp, continuous_3d_dir,\
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
          (python_exe_fp, qiime_dir, mapping_fields, pc_fp, discrete_3d_dir,\
           mapping_fp, params_str)
       
        commands.append([\
          ('Make 3D plots (continuous coloring, %s)' %\
            beta_diversity_metric,continuous_3d_command),\
          ('Make 3D plots (discrete coloring, %s)' %\
            beta_diversity_metric,discrete_3d_command,)])
    
    # Call the command handler on the list of commands
    command_handler(commands,status_update_callback)
    
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
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    qiime_home = qiime_config['qiime_home']
    qiime_dir = qiime_config['qiime_dir']
    
    alpha_diversity_metrics = params['alpha_diversity']['metrics'].split(',')
    
    # Prep the rarefaction command
    min_count, max_count, median_count, mean_count, counts_per_sample =\
     compute_seqs_per_library_stats(open(otu_table_fp,'U'))
    step = int((median_count - min_seqs_per_sample) / num_steps)
    median_count = int(median_count)
    
    rarefaction_dir = '%s/rarefaction/' % output_dir
    try:
        makedirs(rarefaction_dir)
    except OSError:
        pass
    try:
        params_str = get_params_str(params['rarefaction'])
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
         '%s %s/parallel/rarefaction.py -T -i %s -m %s -x %s -s %s -o %s %s' %\
         (python_exe_fp, qiime_dir, otu_table_fp, min_seqs_per_sample, median_count, \
          step, rarefaction_dir, params_str)
    else:
        # Build the rarefaction command
        rarefaction_cmd = \
         '%s %s/rarefaction.py -i %s -m %s -x %s -s %s -o %s %s' %\
         (python_exe_fp, qiime_dir, otu_table_fp, min_seqs_per_sample, median_count, \
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
         "%s %s/parallel/alpha_diversity.py -T -i %s -o %s -t %s %s" %\
         (python_exe_fp, qiime_dir, rarefaction_dir, alpha_diversity_dir, \
          tree_fp, params_str)
    else:  
        # Build the alpha diversity command
        alpha_diversity_cmd = \
         "%s %s/alpha_diversity.py -i %s -o %s -t %s %s" %\
         (python_exe_fp, qiime_dir, rarefaction_dir, alpha_diversity_dir, \
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
     (python_exe_fp, qiime_dir, alpha_diversity_dir, \
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
    for metric in alpha_diversity_metrics:
        input_fp = '%s/%s.txt' % (alpha_collated_dir, metric)
        make_rarefaction_plot_cmd =\
         '%s %s/make_rarefaction_plots.py -m %s -r %s -o %s %s' %\
         (python_exe_fp, qiime_dir, mapping_fp, input_fp, \
          rarefaction_plot_dir, params_str)
        commands.append(\
         [('Rarefaction plot: %s' % metric,make_rarefaction_plot_cmd)])
    
    # Call the command handler on the list of commands
    command_handler(commands,status_update_callback)

    
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
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    qiime_home = qiime_config['qiime_home']
    qiime_dir = qiime_config['qiime_dir']
    
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
     (python_exe_fp, qiime_dir, otu_table_fp, output_dir, params_str)
    commands.append(\
     [('Beta Diversity (%s)' % ', '.join(beta_diversity_metrics), beta_div_cmd)])

    # Prep rarefaction command
    rarefaction_dir = '%s/rarefaction/' % output_dir
    try:
        makedirs(rarefaction_dir)
    except OSError:
        pass
    try:
        params_str = get_params_str(params['rarefaction'])
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
        # Build the parallel rarefaction command
        rarefaction_cmd = \
         '%s %s/parallel/rarefaction.py -T -i %s -m %s -x %s -s 1 -o %s %s' %\
         (python_exe_fp, qiime_dir, otu_table_fp, seqs_per_sample,\
          seqs_per_sample, rarefaction_dir, params_str)
    else:
        # Build the serial rarefaction command
        rarefaction_cmd = \
         '%s %s/rarefaction.py -i %s -m %s -x %s -s 1 -o %s %s' %\
         (python_exe_fp, qiime_dir, otu_table_fp, seqs_per_sample, \
          seqs_per_sample, rarefaction_dir, params_str)
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
            params_str = get_params_str(params['hierarchical_cluster'])
        except KeyError:
            params_str = ''
        # Build the hierarchical clustering command (for full distance matrix)
        hierarchical_cluster_cmd = '%s %s/hierarchical_cluster.py -i %s -o %s %s' %\
         (python_exe_fp, qiime_dir, distance_matrix_fp, master_tree_fp, params_str)
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
             '%s %s/parallel/beta_diversity.py -T -i %s -o %s %s' %\
             (python_exe_fp, qiime_dir, rarefaction_dir, dm_dir, params_str)
        else:
            # Build the serial beta diversity command (for rarefied OTU tables)
            beta_div_rarefied_cmd = \
             '%s %s/beta_diversity.py -i %s -o %s %s' %\
             (python_exe_fp, qiime_dir, rarefaction_dir, dm_dir, params_str)
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
            params_str = get_params_str(params['hierarchical_cluster'])
        except KeyError:
            params_str = ''
        # Build the hierarchical clustering command (for rarefied 
        # distance matrices)
        hierarchical_cluster_cmd =\
         '%s %s/hierarchical_cluster.py -i %s -o %s %s' %\
         (python_exe_fp, qiime_dir, dm_dir, upgma_dir, params_str)
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
         (python_exe_fp, qiime_dir, upgma_dir, master_tree_fp, \
          tree_compare_dir, params_str)
        commands.append(\
         [('Tree compare (%s)' % beta_diversity_metric,\
           tree_compare_cmd)])
           
    # Call the command handler on the list of commands
    command_handler(commands,status_update_callback)
    
    
## End task-specific workflow functions
    