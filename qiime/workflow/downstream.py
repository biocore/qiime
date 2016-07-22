#!/usr/bin/env python
# File created on 20 Feb 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Justin Kuczynski",
               "Jesse Stombaugh", "Yoshiki Vazquez Baeza", "Jai Ram Rideout",
               "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import split, splitext, join
from shutil import rmtree

from biom import load_table
from biom.util import compute_counts_per_sample_stats

from qiime.parse import parse_mapping_file
from qiime.util import create_dir, get_interesting_mapping_fields
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
                                     tree_fp=None,
                                     parallel=False,
                                     logger=None,
                                     suppress_emperor_plots=False,
                                     suppress_md5=False,
                                     status_update_callback=print_to_stdout):
    """ Compute beta diversity distance matrices, run PCoA, and generate emperor plots

        The steps performed by this function are:
         1) Compute a beta diversity distance matrix for each metric
         2) Peform a principal coordinates analysis on the result of step 1
         3) Generate an emperor plot for each result of step 2

    """
    # Prepare some variables for the later steps
    otu_table_dir, otu_table_filename = split(otu_table_fp)
    otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    create_dir(output_dir)
    commands = []
    if logger is None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False

    if not suppress_md5:
        log_input_md5s(logger, [otu_table_fp, mapping_fp, tree_fp])

    mapping_data, mapping_header, mapping_comments =\
        parse_mapping_file(open(mapping_fp, 'U'))

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
            'single_rarefaction.py -i %s -o %s -d %d' %\
            (otu_table_fp, even_sampled_otu_table_fp, sampling_depth)
        commands.append([
            ('Sample OTU table at %d seqs/sample' % sampling_depth,
             single_rarefaction_cmd)])
        otu_table_fp = even_sampled_otu_table_fp
        otu_table_dir, otu_table_filename = split(even_sampled_otu_table_fp)
        otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    try:
        beta_diversity_metrics = params['beta_diversity']['metrics'].split(',')
    except KeyError:
        beta_diversity_metrics = ['weighted_unifrac', 'unweighted_unifrac']

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
            params_str = '%s -t %s ' % (params_str, tree_fp)

        # Build the beta-diversity command
        if parallel:
            # Grab the parallel-specific parameters
            try:
                params_str += get_params_str(params['parallel'])
            except KeyError:
                pass
            beta_div_cmd = 'parallel_beta_diversity.py -i %s -o %s --metrics %s -T %s' %\
                (otu_table_fp, output_dir, beta_diversity_metric, params_str)
            commands.append(
                [('Beta Diversity (%s)' % beta_diversity_metric, beta_div_cmd)])
        else:
            beta_div_cmd = 'beta_diversity.py -i %s -o %s --metrics %s %s' %\
                (otu_table_fp, output_dir, beta_diversity_metric, params_str)
            commands.append(
                [('Beta Diversity (%s)' % beta_diversity_metric, beta_div_cmd)])

        orig_beta_div_fp = '%s/%s_%s.txt' % \
            (output_dir, beta_diversity_metric, otu_table_basename)
        beta_div_fp = '%s/%s_dm.txt' % \
            (output_dir, beta_diversity_metric)
        commands.append(
            [('Rename distance matrix (%s)' % beta_diversity_metric,
              'mv %s %s' % (orig_beta_div_fp, beta_div_fp))])
        dm_fps.append((beta_diversity_metric, beta_div_fp))

        # Prep the principal coordinates command
        pc_fp = '%s/%s_pc.txt' % (output_dir, beta_diversity_metric)
        try:
            params_str = get_params_str(params['principal_coordinates'])
        except KeyError:
            params_str = ''
        # Build the principal coordinates command
        pc_cmd = 'principal_coordinates.py -i %s -o %s %s' %\
            (beta_div_fp, pc_fp, params_str)
        commands.append(
            [('Principal coordinates (%s)' % beta_diversity_metric, pc_cmd)])

        # Generate emperor plots
        if not suppress_emperor_plots:
            # Prep the emperor plots command
            emperor_dir = '%s/%s_emperor_pcoa_plot/' % (output_dir,
                                                        beta_diversity_metric)
            create_dir(emperor_dir)
            try:
                params_str = get_params_str(params['make_emperor'])
            except KeyError:
                params_str = ''
            # Build the continuous-coloring 3d plots command
            emperor_command = \
                'make_emperor.py -i %s -o %s -m %s %s' % (pc_fp,
                                                          emperor_dir,
                                                          mapping_fp,
                                                          params_str)

            commands.append(
                [('Make emperor plots, %s)' % beta_diversity_metric,
                  emperor_command)])

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
                          plot_stderr_and_stddev=False,
                          retain_intermediate_files=True):
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
    if logger is None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False

    if not suppress_md5:
        log_input_md5s(logger, [otu_table_fp, mapping_fp, tree_fp])

    if max_rare_depth is None:
        min_count, max_count, median_count, mean_count, counts_per_sample =\
            compute_counts_per_sample_stats(
                load_table(otu_table_fp))
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
            'parallel_multiple_rarefactions.py -T -i %s -m %s -x %s -s %s -o %s %s' %\
            (otu_table_fp, min_rare_depth, max_rare_depth, step,
             rarefaction_dir, params_str)
    else:
        # Build the rarefaction command
        rarefaction_cmd = \
            'multiple_rarefactions.py -i %s -m %s -x %s -s %s -o %s %s' %\
            (otu_table_fp, min_rare_depth, max_rare_depth, step,
             rarefaction_dir, params_str)
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
            "parallel_alpha_diversity.py -T -i %s -o %s %s" %\
            (rarefaction_dir, alpha_diversity_dir, params_str)
    else:
        # Build the alpha diversity command
        alpha_diversity_cmd = \
            "alpha_diversity.py -i %s -o %s %s" %\
            (rarefaction_dir, alpha_diversity_dir, params_str)

    commands.append(
        [('Alpha diversity on rarefied OTU tables', alpha_diversity_cmd)])

    # Prep the alpha diversity collation command
    alpha_collated_dir = '%s/alpha_div_collated/' % output_dir
    create_dir(alpha_collated_dir)
    try:
        params_str = get_params_str(params['collate_alpha'])
    except KeyError:
        params_str = ''
    # Build the alpha diversity collation command
    alpha_collated_cmd = 'collate_alpha.py -i %s -o %s %s' %\
        (alpha_diversity_dir, alpha_collated_dir, params_str)
    commands.append([('Collate alpha', alpha_collated_cmd)])

    if not retain_intermediate_files:
        commands.append([('Removing intermediate files',
                          'rm -r %s %s' % (rarefaction_dir, alpha_diversity_dir))])
    else:
        commands.append([('Skipping removal of intermediate files.', '')])

    # Prep the make rarefaction plot command(s)
    try:
        params_str = get_params_str(params['make_rarefaction_plots'])
    except KeyError:
        params_str = ''

    if 'std_type' in params['make_rarefaction_plots'] or not plot_stderr_and_stddev:
        rarefaction_plot_dir = '%s/alpha_rarefaction_plots/' % output_dir
        create_dir(rarefaction_plot_dir)

        # Build the make rarefaction plot command(s)
        # for metric in alpha_diversity_metrics:
        make_rarefaction_plot_cmd =\
            'make_rarefaction_plots.py -i %s -m %s -o %s %s' %\
            (alpha_collated_dir, mapping_fp, rarefaction_plot_dir, params_str)
        commands.append(
            [('Rarefaction plot: %s' % 'All metrics', make_rarefaction_plot_cmd)])
    else:
        rarefaction_plot_dir_stddev = '%s/alpha_rarefaction_plots_stddev/' % output_dir
        rarefaction_plot_dir_stderr = '%s/alpha_rarefaction_plots_stderr/' % output_dir
        create_dir(rarefaction_plot_dir_stddev)
        create_dir(rarefaction_plot_dir_stderr)

        # Build the make rarefaction plot command(s)
        # for metric in alpha_diversity_metrics:
        make_rarefaction_plot_cmd =\
            'make_rarefaction_plots.py -i %s -m %s -o %s %s --std_type stddev' %\
            (alpha_collated_dir, mapping_fp, rarefaction_plot_dir_stddev,
             params_str)
        commands.append(
            [('Rarefaction plot: %s' % 'All metrics', make_rarefaction_plot_cmd)])
        make_rarefaction_plot_cmd =\
            'make_rarefaction_plots.py -i %s -m %s -o %s %s --std_type stderr' %\
            (alpha_collated_dir, mapping_fp, rarefaction_plot_dir_stderr,
             params_str)
        commands.append(
            [('Rarefaction plot: %s' % 'All metrics', make_rarefaction_plot_cmd)])

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
    if master_tree is None:
        master_tree = 'full'
    otu_table_dir, otu_table_filename = split(otu_table_fp)
    otu_table_basename, otu_table_ext = splitext(otu_table_filename)
    create_dir(output_dir)
    commands = []
    if logger is None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False

    if not suppress_md5:
        log_input_md5s(logger, [otu_table_fp, mapping_fp, tree_fp])

    try:
        beta_diversity_metrics = params['beta_diversity']['metrics'].split(',')
    except KeyError:
        beta_diversity_metrics = ['weighted_unifrac', 'unweighted_unifrac']

    # Prep the beta-diversity command
    full_table_bdiv_dir = join(output_dir, 'unrarefied_bdiv')
    try:
        params_str = get_params_str(params['beta_diversity'])
    except KeyError:
        params_str = ''
    if tree_fp:
        params_str = '%s -t %s' % (params_str, tree_fp)
    # Build the beta-diversity command
    beta_div_cmd = 'beta_diversity.py -i %s -o %s %s' %\
        (otu_table_fp, full_table_bdiv_dir, params_str)
    commands.append(
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
        'multiple_rarefactions_even_depth.py -i %s -d %d -o %s %s' %\
        (otu_table_fp, seqs_per_sample, rarefaction_dir, params_str)
    commands.append([('Rarefaction', rarefaction_cmd)])

    # Begin iterating over beta diversity distance metrics, if more than one
    # was provided
    for beta_diversity_metric in beta_diversity_metrics:
        metric_output_dir = '%s/%s/' % (output_dir, beta_diversity_metric)
        distance_matrix_fp = '%s/%s_%s.txt' % \
            (full_table_bdiv_dir, beta_diversity_metric, otu_table_basename)

        # Prep the hierarchical clustering command (for full distance matrix)
        full_tree_fp = '%s/%s_%s_upgma.tre' % (full_table_bdiv_dir,
                                               otu_table_basename,
                                               beta_diversity_metric)
        try:
            params_str = get_params_str(params['upgma_cluster'])
        except KeyError:
            params_str = ''
        # Build the hierarchical clustering command (for full distance matrix)
        hierarchical_cluster_cmd = 'upgma_cluster.py -i %s -o %s %s' %\
            (distance_matrix_fp, full_tree_fp, params_str)
        commands.append(
            [('UPGMA on full distance matrix: %s' % beta_diversity_metric,
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
            params_str = '%s -t %s' % (params_str, tree_fp)
        if parallel:
            params_str += ' %s' % get_params_str(params['parallel'])
            # Build the parallel beta diversity command (for rarefied OTU
            # tables)
            beta_div_rarefied_cmd = \
                'parallel_beta_diversity.py -T -i %s -o %s %s' %\
                (rarefaction_dir, dm_dir, params_str)
        else:
            # Build the serial beta diversity command (for rarefied OTU tables)
            beta_div_rarefied_cmd = \
                'beta_diversity.py -i %s -o %s %s' %\
                (rarefaction_dir, dm_dir, params_str)
        commands.append(
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
            'upgma_cluster.py -i %s -o %s %s' % (dm_dir, upgma_dir, params_str)
        commands.append(
            [('UPGMA on rarefied distance matrix (%s)' % beta_diversity_metric,
              hierarchical_cluster_cmd)])

        # Build the consensus tree command
        consensus_tree_cmd =\
            'consensus_tree.py -i %s -o %s %s' %\
            (upgma_dir, metric_output_dir + "/rare_upgma_consensus.tre",
             params_str)
        commands.append(
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
            master_tree_fp = metric_output_dir + "/rare_upgma_consensus.tre"
        else:
            raise RuntimeError(
                'master tree method "%s" not found' %
                (master_tree,))
        tree_compare_cmd = 'tree_compare.py -s %s -m %s -o %s %s' %\
            (upgma_dir, master_tree_fp, tree_compare_dir, params_str)
        commands.append(
            [('Tree compare (%s)' % beta_diversity_metric, tree_compare_cmd)])

        # Prep the PCoA command
        pcoa_dir = '%s/pcoa/' % metric_output_dir
        create_dir(pcoa_dir)
        try:
            params_str = get_params_str(params['principal_coordinates'])
        except KeyError:
            params_str = ''
        # Build the PCoA command
        pcoa_cmd = 'principal_coordinates.py -i %s -o %s %s' %\
            (dm_dir, pcoa_dir, params_str)
        commands.append(
            [('Principal coordinates (%s)' % beta_diversity_metric, pcoa_cmd)])

        # Prep the emperor plots command
        emperor_dir = '%s/emperor_pcoa_plots/' % metric_output_dir
        create_dir(emperor_dir)
        try:
            params_str = get_params_str(params['make_emperor'])
        except KeyError:
            params_str = ''
        emperor_cmd = 'make_emperor.py -i %s -o %s -m %s %s' %\
            (pcoa_dir, emperor_dir, mapping_fp, params_str)
        commands.append(
            [('emperor plots (%s)' % beta_diversity_metric, emperor_cmd)])

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

    if logger is None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False

    if not suppress_md5:
        log_input_md5s(logger, [otu_table_fp, mapping_fp])

    # if mapping category not passed via command-line,
    # check if it is passed in params file
    if not mapping_cat:
        try:
            mapping_cat = params['collapse_samples']['collapse_fields']
        except:
            mapping_cat = None

    try:
        params_str = get_params_str(params['collapse_samples'])
        # Need to remove the mapping category option, since it is defined above.
        # Using this method since we don't want to change the params dict
        split_params = params_str.split('--')
        updated_params_str = []
        for i in split_params:
            if not i.startswith('collapse_fields'):
                updated_params_str.append(i)
        params_str = '--'.join(updated_params_str)
    except:
        params_str = ''

    if mapping_cat:
        base_filename = mapping_cat.replace(' ', '-').replace(',','')
        output_biom_fp = join(
            output_dir, '%s_otu_table.biom' % base_filename)
        output_map_fp = join(
            output_dir, '%s_map.txt' % base_filename)
        # Build the collapse samples command
        collapse_samples_cmd = \
            "collapse_samples.py -m %s -b %s --output_biom_fp %s --output_mapping_fp %s --collapse_fields '%s' %s" %\
            (mapping_fp, otu_table_fp, output_biom_fp, output_map_fp, mapping_cat, params_str)

        commands.append(
            [('Collapse samples in OTU table by categories', collapse_samples_cmd)])

        otu_table_fp = output_biom_fp

    # Build the sort OTU table command
    if sort:
        # Prep the sort_otu_table command
        try:
            params_str = get_params_str(params['sort_otu_table'])
        except:
            params_str = ''

        # define output otu table
        sorted_fp = join(output_dir,
                         splitext(split(otu_table_fp)[-1])[0] + '_sorted.biom')

        if mapping_cat or params_str == '':
            # for this case we don't have a collapsed mapping file so must
            # handle separately
            sort_otu_table_cmd = \
                "sort_otu_table.py -i %s -o %s" % (otu_table_fp, sorted_fp)
        else:
            sort_otu_table_cmd = \
                "sort_otu_table.py -i %s -o %s -m %s %s" %\
                (otu_table_fp, sorted_fp, mapping_fp, params_str)

        commands.append([('Sort OTU Table', sort_otu_table_cmd)])

        # redefine otu_table_fp to use
        otu_table_fp = sorted_fp

    # Prep the summarize taxonomy command
    try:
        params_str = get_params_str(params['summarize_taxa'])
    except:
        params_str = ''

    try:
        sum_taxa_levels = params['summarize_taxa']['level']
    except:
        sum_taxa_levels = None

    # Build the summarize taxonomy command
    summarize_taxa_cmd = 'summarize_taxa.py -i %s -o %s %s' %\
        (otu_table_fp, output_dir, params_str)

    commands.append([('Summarize Taxonomy', summarize_taxa_cmd)])

    sum_taxa_fps = []

    if sum_taxa_levels:
        basename = join(output_dir, splitext(split(otu_table_fp)[-1])[0])
        for i in sum_taxa_levels.split(','):
            sum_taxa_fps.append(basename + '_L%s.txt' % (str(i)))
    else:
        basename = join(output_dir, splitext(split(otu_table_fp)[-1])[0])
        # this is the default levels from summarize_taxa, but cannot import
        # script to get these values
        for i in [2, 3, 4, 5, 6]:
            sum_taxa_fps.append(basename + '_L%s.txt' % (str(i)))

    # Prep the plot taxa summary plot command(s)
    taxa_summary_plots_dir = '%s/taxa_summary_plots/' % output_dir
    create_dir(taxa_summary_plots_dir)

    try:
        params_str = get_params_str(params['plot_taxa_summary'])
    except:
        params_str = ''
    # Build the plot taxa summary plot command(s)

    plot_taxa_summary_cmd =\
        'plot_taxa_summary.py -i %s -o %s %s' %\
        (','.join(sum_taxa_fps), taxa_summary_plots_dir, params_str)

    commands.append(
        [('Plot Taxonomy Summary', plot_taxa_summary_cmd)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)
