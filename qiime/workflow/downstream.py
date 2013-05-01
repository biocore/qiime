#!/usr/bin/env python
# File created on 20 Feb 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso",
               "Kyle Bittinger",
               "Justin Kuczynski",
               "Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os.path import split, splitext, join
from biom.parse import parse_biom_table
from biom.util import compute_counts_per_sample_stats
from qiime.parse import parse_mapping_file
from qiime.util import (get_qiime_scripts_dir,
                        create_dir,
                        get_interesting_mapping_fields)
from qiime.workflow.util import (print_to_stdout,
                                 generate_log_fp,
                                 WorkflowLogger,
                                 log_input_md5s,
                                 get_params_str)


def run_beta_diversity_through_plots(otu_table_fp, 
                                     mapping_fp,
                                     output_dir,
                                     command_handler,
                                     params,
                                     qiime_config,
                                     color_by_interesting_fields_only=True,
                                     sampling_depth=None,
                                     histogram_categories=None,
                                     tree_fp=None,
                                     parallel=False,
                                     logger=None,
                                     suppress_3d_plots=False,
                                     suppress_2d_plots=False,
                                     suppress_md5=False,
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
    
    if not suppress_md5:
        log_input_md5s(logger,[otu_table_fp,mapping_fp,tree_fp])
    
    mapping_data, mapping_header, mapping_comments =\
     parse_mapping_file(open(mapping_fp,'U'))
    if histogram_categories:
        invalid_categories = set(histogram_categories) - set(mapping_header)
        if invalid_categories:
            raise ValueError,\
             "Invalid histogram categories - these must exactly match "+\
             "mapping file column headers: %s" % (' '.join(invalid_categories))
    # Get the interesting mapping fields to color by -- if none are
    # interesting, take all of them. Interesting is defined as those
    # which have greater than one value and fewer values than the number 
    # of samples
    if color_by_interesting_fields_only:
        mapping_fields =\
          get_interesting_mapping_fields(mapping_data, mapping_header) or\
          mapping_header
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
        
        
        orig_beta_div_fp = '%s/%s_%s.txt' % \
         (output_dir, beta_diversity_metric, otu_table_basename)
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
            create_dir(continuous_3d_dir)
            try:
                params_str = get_params_str(params['make_3d_plots'])
            except KeyError:
                params_str = ''
            # Build the continuous-coloring 3d plots command
            continuous_3d_command = \
             '%s %s/make_3d_plots.py -p %s -i %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, prefs_fp, pc_fp, continuous_3d_dir,
               mapping_fp, params_str)
    
            # Prep the discrete-coloring 3d plots command
            discrete_3d_dir = '%s/%s_3d_discrete/' %\
             (output_dir, beta_diversity_metric)
            create_dir(discrete_3d_dir)
            try:
                params_str = get_params_str(params['make_3d_plots'])
            except KeyError:
                params_str = ''
            # Build the discrete-coloring 3d plots command
            discrete_3d_command = \
             '%s %s/make_3d_plots.py -b "%s" -i %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, mapping_fields, pc_fp, discrete_3d_dir,
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
            create_dir(continuous_2d_dir)
            try:
                params_str = get_params_str(params['make_2d_plots'])
            except KeyError:
                params_str = ''
            # Build the continuous-coloring 3d plots command
            continuous_2d_command = \
             '%s %s/make_2d_plots.py -p %s -i %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, prefs_fp, pc_fp, continuous_2d_dir,
               mapping_fp, params_str)
               
            # Prep the discrete-coloring 3d plots command
            discrete_2d_dir = '%s/%s_2d_discrete/' %\
             (output_dir, beta_diversity_metric)
            create_dir(discrete_2d_dir)
            try:
                params_str = get_params_str(params['make_2d_plots'])
            except KeyError:
                params_str = ''
            # Build the discrete-coloring 2d plots command
            discrete_2d_command = \
             '%s %s/make_2d_plots.py -b "%s" -i %s -o %s -m %s %s' %\
              (python_exe_fp, script_dir, mapping_fields, pc_fp, discrete_2d_dir,
               mapping_fp, params_str)
       
            commands.append([\
              ('Make 2D plots (continuous coloring, %s)' %\
                beta_diversity_metric,continuous_2d_command),\
              ('Make 2D plots (discrete coloring, %s)' %\
                beta_diversity_metric,discrete_2d_command,)])
                
        if histogram_categories:
            # Prep the discrete-coloring 3d plots command
            histograms_dir = '%s/%s_histograms/' %\
             (output_dir, beta_diversity_metric)
            create_dir(histograms_dir)
            try:
                params_str = get_params_str(params['make_distance_histograms'])
            except KeyError:
                params_str = ''
            # Build the make_distance_histograms command
            distance_histograms_command = \
             '%s %s/make_distance_histograms.py -d %s -o %s -m %s -f "%s" %s' %\
              (python_exe_fp, script_dir, beta_div_fp, 
               histograms_dir, mapping_fp, 
               ','.join(histogram_categories), params_str)
       
            commands.append([\
              ('Make Distance Histograms (%s)' %\
                beta_diversity_metric,distance_histograms_command)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)
    
    return dm_fps

def run_alpha_rarefaction(otu_table_fp, 
                          mapping_fp,
                          output_dir,
                          command_handler,
                          params,
                          qiime_config,
                          tree_fp=None,
                          num_steps=10,
                          parallel=False,
                          logger=None,
                          min_rare_depth=10,
                          max_rare_depth=None,
                          suppress_md5=False,
                          status_update_callback=print_to_stdout,
                          plot_stderr_and_stddev=False):
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

    if not suppress_md5:
        log_input_md5s(logger,[otu_table_fp,mapping_fp,tree_fp])
    
    if max_rare_depth == None:
        min_count, max_count, median_count, mean_count, counts_per_sample =\
         compute_counts_per_sample_stats(parse_biom_table(open(otu_table_fp,'U')))
        max_rare_depth = median_count
    step = int((max_rare_depth - min_rare_depth) / num_steps) or 1
    max_rare_depth = int(max_rare_depth)
    
    rarefaction_dir = '%s/rarefaction/' % output_dir
    create_dir(rarefaction_dir)
    try:
        params_str = get_params_str(params['multiple_rarefactions'])
    except KeyError:
        params_str = ''
    if parallel:
        params_str += ' %s' % get_params_str(params['parallel'])        
        # Build the rarefaction command
        rarefaction_cmd = \
         '%s %s/parallel_multiple_rarefactions.py -T -i %s -m %s -x %s -s %s -o %s %s' %\
         (python_exe_fp, script_dir, otu_table_fp, min_rare_depth, max_rare_depth,
          step, rarefaction_dir, params_str)
    else:
        # Build the rarefaction command
        rarefaction_cmd = \
         '%s %s/multiple_rarefactions.py -i %s -m %s -x %s -s %s -o %s %s' %\
         (python_exe_fp, script_dir, otu_table_fp, min_rare_depth, max_rare_depth,
          step, rarefaction_dir, params_str)
    commands.append([('Alpha rarefaction', rarefaction_cmd)])
    
    # Prep the alpha diversity command
    alpha_diversity_dir = '%s/alpha_div/' % output_dir
    create_dir(alpha_diversity_dir)
    try:
        params_str = get_params_str(params['alpha_diversity'])
    except KeyError:
        params_str = ''
    if tree_fp:
        params_str += ' -t %s' % tree_fp
    if parallel:
        params_str += ' %s' % get_params_str(params['parallel'])   
        # Build the alpha diversity command
        alpha_diversity_cmd = \
         "%s %s/parallel_alpha_diversity.py -T -i %s -o %s %s" %\
         (python_exe_fp, script_dir, rarefaction_dir, alpha_diversity_dir,
          params_str)
    else:  
        # Build the alpha diversity command
        alpha_diversity_cmd = \
         "%s %s/alpha_diversity.py -i %s -o %s %s" %\
         (python_exe_fp, script_dir, rarefaction_dir, alpha_diversity_dir,
          params_str)

    commands.append(\
     [('Alpha diversity on rarefied OTU tables',alpha_diversity_cmd)])
     
    # Prep the alpha diversity collation command
    alpha_collated_dir = '%s/alpha_div_collated/' % output_dir
    create_dir(alpha_collated_dir)
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
    try:
        params_str = get_params_str(params['make_rarefaction_plots'])
    except KeyError:
        params_str = ''
    
    if 'std_type' in params['make_rarefaction_plots'] or not plot_stderr_and_stddev:
        rarefaction_plot_dir = '%s/alpha_rarefaction_plots/' % output_dir
        create_dir(rarefaction_plot_dir)
        
        # Build the make rarefaction plot command(s)
        #for metric in alpha_diversity_metrics:
        make_rarefaction_plot_cmd =\
             '%s %s/make_rarefaction_plots.py -i %s -m %s -o %s %s' %\
             (python_exe_fp, script_dir, alpha_collated_dir, mapping_fp,
              rarefaction_plot_dir, params_str)
        commands.append(\
             [('Rarefaction plot: %s' % 'All metrics',make_rarefaction_plot_cmd)])
    else:
        rarefaction_plot_dir_stddev = '%s/alpha_rarefaction_plots_stddev/' % output_dir
        rarefaction_plot_dir_stderr = '%s/alpha_rarefaction_plots_stderr/' % output_dir
        create_dir(rarefaction_plot_dir_stddev)
        create_dir(rarefaction_plot_dir_stderr)
        
        # Build the make rarefaction plot command(s)
        # for metric in alpha_diversity_metrics:
        make_rarefaction_plot_cmd =\
             '%s %s/make_rarefaction_plots.py -i %s -m %s -o %s %s --std_type stddev' %\
             (python_exe_fp, script_dir, alpha_collated_dir, mapping_fp,
              rarefaction_plot_dir_stddev, params_str)
        commands.append(\
             [('Rarefaction plot: %s' % 'All metrics',make_rarefaction_plot_cmd)])
        make_rarefaction_plot_cmd =\
             '%s %s/make_rarefaction_plots.py -i %s -m %s -o %s %s --std_type stderr' %\
             (python_exe_fp, script_dir, alpha_collated_dir, mapping_fp,
              rarefaction_plot_dir_stderr, params_str)
        commands.append(\
             [('Rarefaction plot: %s' % 'All metrics',make_rarefaction_plot_cmd)])
   
    
    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)


run_qiime_alpha_rarefaction = run_alpha_rarefaction

def run_jackknifed_beta_diversity(otu_table_fp,
                                  tree_fp,
                                  seqs_per_sample,
                                  output_dir,
                                  command_handler,
                                  params,
                                  qiime_config,
                                  mapping_fp,
                                  parallel=False,
                                  logger=None,
                                  suppress_md5=False,
                                  status_update_callback=print_to_stdout,
                                  master_tree=None):
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
    
    if not suppress_md5:
        log_input_md5s(logger,[otu_table_fp,mapping_fp,tree_fp])
    
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
    create_dir(rarefaction_dir)
    try:
        params_str = get_params_str(params['multiple_rarefactions_even_depth'])
    except KeyError:
        params_str = ''
    # Build the rarefaction command
    rarefaction_cmd = \
     '%s %s/multiple_rarefactions_even_depth.py -i %s -d %d -o %s %s' %\
     (python_exe_fp, script_dir, otu_table_fp, seqs_per_sample,
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
        create_dir(dm_dir)
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
         [('Beta diversity on rarefied OTU tables (%s)' % beta_diversity_metric,
           beta_div_rarefied_cmd)])

        # Prep the hierarchical clustering command (for rarefied 
        # distance matrices)
        upgma_dir = '%s/rare_upgma/' % metric_output_dir
        create_dir(upgma_dir)

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
         [('UPGMA on rarefied distance matrix (%s)' % beta_diversity_metric,
           hierarchical_cluster_cmd)])
        

        # Build the consensus tree command
        consensus_tree_cmd =\
         '%s %s/consensus_tree.py -i %s -o %s %s' %\
         (python_exe_fp, script_dir, upgma_dir, upgma_dir + "/consensus.tre",
            params_str)
        commands.append(\
         [('consensus on rarefied distance matrices (%s)' % beta_diversity_metric,
           consensus_tree_cmd)])
           
           
        # Prep the tree compare command
        tree_compare_dir = '%s/upgma_cmp/' % metric_output_dir
        create_dir(tree_compare_dir)
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
         (python_exe_fp, script_dir, upgma_dir, master_tree_fp, 
          tree_compare_dir, params_str)
        commands.append(\
         [('Tree compare (%s)' % beta_diversity_metric, tree_compare_cmd)])
           
        # Prep the PCoA command
        pcoa_dir = '%s/pcoa/' % metric_output_dir
        create_dir(pcoa_dir)
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
        create_dir(plots_2d_dir)
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
        create_dir(plots_3d_dir)
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
    

def run_summarize_taxa_through_plots(otu_table_fp, 
                                     mapping_fp,
                                     output_dir,
                                     mapping_cat,
                                     sort,
                                     command_handler,
                                     params,
                                     qiime_config,
                                     logger=None, 
                                     suppress_md5=False,
                                     status_update_callback=print_to_stdout):
    """ Run the data preparation for summarizing taxonomies and generating plots
    
        The steps performed by this function are:
          1) Summarize OTU by Category
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
    
    if not suppress_md5:
        log_input_md5s(logger,[otu_table_fp,mapping_fp])
    
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
        output_fp=join(output_dir,'%s_otu_table.biom' % (mapping_cat.replace(' ','-')))
        # Build the summarize otu by category command
        summarize_otu_by_cat_cmd = \
         "%s %s/summarize_otu_by_cat.py -i %s -c %s -o %s -m '%s' %s" %\
         (python_exe_fp, script_dir, mapping_fp, otu_table_fp, output_fp,
          mapping_cat, params_str)
        
        commands.append(\
         [('Summarize OTU table by Category',summarize_otu_by_cat_cmd)])
         
        otu_table_fp=output_fp
    
    # Build the sort OTU table command
    if sort:
        # Prep the sort_otu_table command
        try:
            params_str = get_params_str(params['sort_otu_table'])
        except:
            params_str = ''
            
        # define output otu table
        sorted_fp=join(output_dir,
                       splitext(split(otu_table_fp)[-1])[0]+'_sorted.biom')
        
        if mapping_cat or params_str=='':
            # for this case we don't have a collapsed mapping file so must
            # handle separately
            sort_otu_table_cmd = \
             "%s %s/sort_otu_table.py -i %s -o %s" %\
             (python_exe_fp, script_dir, otu_table_fp, sorted_fp)
        else:
            sort_otu_table_cmd = \
             "%s %s/sort_otu_table.py -i %s -o %s -m %s %s" %\
             (python_exe_fp, script_dir, otu_table_fp, sorted_fp,
              mapping_fp, params_str)
        
        commands.append([('Sort OTU Table',sort_otu_table_cmd)])

        # redefine otu_table_fp to use
        otu_table_fp=sorted_fp
    
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
     (python_exe_fp, script_dir, otu_table_fp, output_dir, params_str)
    
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
    create_dir(taxa_summary_plots_dir)
        
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
