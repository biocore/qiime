#!/usr/bin/env python
# File created on 09 Feb 2013
from __future__ import division
import re
from glob import glob
from os.path import split, splitext, exists, join
from qiime.parse import (parse_qiime_parameters,
                         parse_mapping_file_to_dict)
from qiime.util import (create_dir,
                        MetadataMap)
from qiime.workflow.downstream import (
    run_beta_diversity_through_plots,
    run_alpha_rarefaction,
    run_summarize_taxa_through_plots)
from qiime.workflow.util import (print_to_stdout,
                                 generate_log_fp,
                                 WorkflowLogger,
                                 log_input_md5s,
                                 call_commands_serially,
                                 get_params_str)

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

#_index_headers includes key/value pairs use in pick_open_reference_otus.py
_index_headers = {
    "run_summary": "Run summary data",
    "beta_diversity_even": "Beta diversity results (even sampling: %d)",
    "alpha_diversity": "Alpha diversity results",
    "taxa_summary": "Taxonomic summary results",
    "taxa_summary_categorical": "Taxonomic summary results (by %s)",
    "group_significance": "Group significance results",
    "failures_fp": "Failures",
    "otu_maps": "OTU maps",
    "sequences": "Sequences",
    "otu_tables": "OTU tables",
    "trees": "Trees",
    "taxa_assignments": "Taxonomy assignments"}


def format_index_link(link_description, relative_path):
    return '<td>%s</td><td> <a href="%s" target="_blank">%s</a></td>' % (link_description,
                                                                         re.sub(
                                                                             '/+',
                                                                             '/',
                                                                             relative_path),
                                                                         split(relative_path)[1])


def generate_index_page(index_links,
                        index_fp,
                        order=[_index_headers['run_summary']]):
    """ generate the top-level index page """
    # get containing directory for index_fp
    top_level_dir = split(split(index_fp)[0])[1]
    index_page_header = get_index_page_header()
    index_lines = [index_page_header]
    d = {}
    for e in index_links:
        try:
            d[e[2]].append((e[0], e[1]))
        except KeyError:
            d[e[2]] = [(e[0], e[1])]
    index_lines.append('<table border=1>\n')

    # Determine the order the data should be presented in. This should be
    # the order that the user requested, followed by any categories that
    # the user didn't include in the order parameter.
    ordered_table_entries = order + [k for k in d if k not in order]
    for k in ordered_table_entries:
        v = d[k]
        index_lines.append(
            '<tr colspan=2 align=center bgcolor=#e8e8e8><td colspan=2 align=center>%s</td></tr>\n' % k)
        for description, path in v:
            # if path starts with top_level_dir, replace it
            # with ./
            path = re.sub('^.*%s\/' % top_level_dir, './', path)
            index_lines.append(
                '<tr>%s</tr>\n' %
                format_index_link(
                    description,
                    path))
    index_lines.append('</table>\n')

    index_page_footer = get_index_page_footer()
    index_lines.append(index_page_footer)

    open(index_fp, 'w').write(''.join(index_lines))


def get_index_page_header():
    return """<html>
<head><title>QIIME results</title></head>
<body>
<a href="http://www.qiime.org" target="_blank"><img src=\"http://qiime.org/_static/wordpressheader.png\" alt="www.qiime.org""/></a><p>
"""


def get_index_page_footer():
    return """<p><b>Need help?</b> See <a href="http://help.qiime.org" target="_blank">http://help.qiime.org</a>.
</body></html>"""


def run_core_diversity_analyses(
        biom_fp,
        mapping_fp,
        sampling_depth,
        output_dir,
        qiime_config,
        command_handler=call_commands_serially,
        tree_fp=None,
        params=None,
        categories=None,
        arare_min_rare_depth=10,
        arare_num_steps=10,
        parallel=False,
        suppress_taxa_summary=False,
        suppress_beta_diversity=False,
        suppress_alpha_diversity=False,
        suppress_group_significance=False,
        status_update_callback=print_to_stdout):
    """
    """
    if categories is not None:
        # Validate categories provided by the users
        mapping_data, mapping_comments = \
            parse_mapping_file_to_dict(open(mapping_fp, 'U'))
        metadata_map = MetadataMap(mapping_data, mapping_comments)
        for c in categories:
            if c not in metadata_map.CategoryNames:
                raise ValueError("Category '%s' is not a column header "
                                 "in your mapping file. "
                                 "Categories are case and white space sensitive. Valid "
                                 "choices are: (%s)" % (c, ', '.join(metadata_map.CategoryNames)))
            if metadata_map.hasSingleCategoryValue(c):
                raise ValueError("Category '%s' contains only one value. "
                                 "Categories analyzed here require at least two values." % c)

    else:
        categories = []
    comma_separated_categories = ','.join(categories)
    # prep some variables
    if params is None:
        params = parse_qiime_parameters([])

    create_dir(output_dir)
    index_fp = '%s/index.html' % output_dir
    index_links = []
    commands = []

    # begin logging
    old_log_fps = glob(join(output_dir, 'log_20*txt'))
    log_fp = generate_log_fp(output_dir)
    index_links.append(
        ('Master run log',
         log_fp,
         _index_headers['run_summary']))
    for old_log_fp in old_log_fps:
        index_links.append(
            ('Previous run log',
             old_log_fp,
             _index_headers['run_summary']))
    logger = WorkflowLogger(log_fp,
                            params=params,
                            qiime_config=qiime_config)
    input_fps = [biom_fp, mapping_fp]
    if tree_fp is not None:
        input_fps.append(tree_fp)
    log_input_md5s(logger, input_fps)

    # run 'biom summarize-table' on input BIOM table
    try:
        params_str = get_params_str(params['biom-summarize-table'])
    except KeyError:
        params_str = ''
    biom_table_stats_output_fp = '%s/biom_table_summary.txt' % output_dir
    if not exists(biom_table_stats_output_fp):
        biom_table_summary_cmd = \
            "biom summarize-table -i %s -o %s %s" % \
            (biom_fp, biom_table_stats_output_fp, params_str)
        commands.append([('Generate BIOM table summary',
                          biom_table_summary_cmd)])
    else:
        logger.write("Skipping 'biom summarize-table' as %s exists.\n\n"
                     % biom_table_stats_output_fp)
    index_links.append(('BIOM table statistics',
                        biom_table_stats_output_fp,
                        _index_headers['run_summary']))

    # filter samples with fewer observations than the requested sampling_depth.
    # since these get filtered for some analyses (eg beta diversity after
    # even sampling) it's useful to filter them here so they're filtered
    # from all analyses.
    filtered_biom_fp = "%s/table_mc%d.biom" % (output_dir, sampling_depth)
    if not exists(filtered_biom_fp):
        filter_samples_cmd = "filter_samples_from_otu_table.py -i %s -o %s -n %d" %\
            (biom_fp, filtered_biom_fp, sampling_depth)
        commands.append(
            [('Filter low sequence count samples from table (minimum sequence count: %d)' % sampling_depth,
              filter_samples_cmd)])
    else:
        logger.write("Skipping filter_samples_from_otu_table.py as %s exists.\n\n"
                     % filtered_biom_fp)
    biom_fp = filtered_biom_fp

    # rarify the BIOM table to sampling_depth
    rarefied_biom_fp = "%s/table_even%d.biom" % (output_dir, sampling_depth)
    if not exists(rarefied_biom_fp):
        single_rarefaction_cmd = "single_rarefaction.py -i %s -o %s -d %d" %\
            (biom_fp, rarefied_biom_fp, sampling_depth)
        commands.append(
            [('Rarify the OTU table to %d sequences/sample' % sampling_depth,
              single_rarefaction_cmd)])
    else:
        logger.write("Skipping single_rarefaction.py as %s exists.\n\n"
                     % rarefied_biom_fp)

    # run initial commands and reset the command list
    if len(commands) > 0:
        command_handler(commands,
                        status_update_callback,
                        logger,
                        close_logger_on_success=False)
        commands = []

    if not suppress_beta_diversity:
        bdiv_even_output_dir = '%s/bdiv_even%d/' % (output_dir, sampling_depth)
        # Need to check for the existence of any distance matrices, since the user
        # can select which will be generated.
        existing_dm_fps = glob('%s/*_dm.txt' % bdiv_even_output_dir)
        if len(existing_dm_fps) == 0:
            even_dm_fps = run_beta_diversity_through_plots(
                otu_table_fp=rarefied_biom_fp,
                mapping_fp=mapping_fp,
                output_dir=bdiv_even_output_dir,
                command_handler=command_handler,
                params=params,
                qiime_config=qiime_config,
                # Note: we pass sampling depth=None here as
                # we rarify the BIOM table above and pass that
                # in here.
                sampling_depth=None,
                tree_fp=tree_fp,
                parallel=parallel,
                logger=logger,
                suppress_md5=True,
                status_update_callback=status_update_callback)
        else:
            logger.write("Skipping beta_diversity_through_plots.py as %s exist(s).\n\n"
                         % ', '.join(existing_dm_fps))
            even_dm_fps = [(split(fp)[1].strip('_dm.txt'), fp)
                           for fp in existing_dm_fps]

        # Get make_distance_boxplots parameters
        try:
            params_str = get_params_str(params['make_distance_boxplots'])
        except KeyError:
            params_str = ''

        for bdiv_metric, dm_fp in even_dm_fps:
            for category in categories:
                boxplots_output_dir = '%s/%s_boxplots/' % (bdiv_even_output_dir,
                                                           bdiv_metric)
                plot_output_fp = '%s/%s_Distances.pdf' % (boxplots_output_dir,
                                                          category)
                stats_output_fp = '%s/%s_Stats.txt' % (boxplots_output_dir,
                                                       category)
                if not exists(plot_output_fp):
                    boxplots_cmd = \
                        'make_distance_boxplots.py -d %s -f %s -o %s -m %s -n 999 %s' %\
                        (dm_fp, category, boxplots_output_dir,
                         mapping_fp, params_str)
                    commands.append([('Boxplots (%s)' % category,
                                      boxplots_cmd)])
                else:
                    logger.write("Skipping make_distance_boxplots.py for %s as %s exists.\n\n"
                                 % (category, plot_output_fp))
                index_links.append(('Distance boxplots (%s)' % bdiv_metric,
                                    plot_output_fp,
                                    _index_headers['beta_diversity_even'] % sampling_depth))
                index_links.append(
                    ('Distance boxplots statistics (%s)' % bdiv_metric,
                     stats_output_fp,
                     _index_headers['beta_diversity_even'] % sampling_depth))

            index_links.append(('PCoA plot (%s)' % bdiv_metric,
                                '%s/%s_emperor_pcoa_plot/index.html' %
                                (bdiv_even_output_dir, bdiv_metric),
                                _index_headers['beta_diversity_even'] % sampling_depth))
            index_links.append(('Distance matrix (%s)' % bdiv_metric,
                                '%s/%s_dm.txt' %
                                (bdiv_even_output_dir, bdiv_metric),
                                _index_headers['beta_diversity_even'] % sampling_depth))
            index_links.append(
                ('Principal coordinate matrix (%s)' % bdiv_metric,
                 '%s/%s_pc.txt' %
                 (bdiv_even_output_dir, bdiv_metric),
                 _index_headers['beta_diversity_even'] % sampling_depth))

    if not suppress_alpha_diversity:
        # Alpha rarefaction workflow
        arare_full_output_dir = '%s/arare_max%d/' % (output_dir,
                                                     sampling_depth)
        rarefaction_plots_output_fp = \
            '%s/alpha_rarefaction_plots/rarefaction_plots.html' % arare_full_output_dir
        if not exists(rarefaction_plots_output_fp):
            run_alpha_rarefaction(
                otu_table_fp=biom_fp,
                mapping_fp=mapping_fp,
                output_dir=arare_full_output_dir,
                command_handler=command_handler,
                params=params,
                qiime_config=qiime_config,
                tree_fp=tree_fp,
                num_steps=arare_num_steps,
                parallel=parallel,
                logger=logger,
                min_rare_depth=arare_min_rare_depth,
                max_rare_depth=sampling_depth,
                suppress_md5=True,
                status_update_callback=status_update_callback,
                retain_intermediate_files=False)
        else:
            logger.write("Skipping alpha_rarefaction.py as %s exists.\n\n"
                         % rarefaction_plots_output_fp)

        index_links.append(('Alpha rarefaction plots',
                            rarefaction_plots_output_fp,
                            _index_headers['alpha_diversity']))

        collated_alpha_diversity_fps = \
            glob('%s/alpha_div_collated/*txt' % arare_full_output_dir)
        try:
            params_str = get_params_str(params['compare_alpha_diversity'])
        except KeyError:
            params_str = ''

        if len(categories) > 0:
            for collated_alpha_diversity_fp in collated_alpha_diversity_fps:
                alpha_metric = splitext(
                    split(collated_alpha_diversity_fp)[1])[0]
                compare_alpha_output_dir = '%s/compare_%s' % \
                    (arare_full_output_dir, alpha_metric)
                if not exists(compare_alpha_output_dir):
                    compare_alpha_cmd = \
                        'compare_alpha_diversity.py -i %s -m %s -c %s -o %s -n 999 %s' %\
                        (collated_alpha_diversity_fp,
                         mapping_fp,
                         comma_separated_categories,
                         compare_alpha_output_dir,
                         params_str)
                    commands.append(
                        [('Compare alpha diversity (%s)' % alpha_metric,
                          compare_alpha_cmd)])
                    for category in categories:
                        alpha_comparison_stat_fp = '%s/%s_stats.txt' % \
                            (compare_alpha_output_dir, category)
                        alpha_comparison_boxplot_fp = '%s/%s_boxplots.pdf' % \
                            (compare_alpha_output_dir, category)
                        index_links.append(
                            ('Alpha diversity statistics (%s, %s)' % (category, alpha_metric),
                             alpha_comparison_stat_fp,
                             _index_headers['alpha_diversity']))
                        index_links.append(
                            ('Alpha diversity boxplots (%s, %s)' % (category, alpha_metric),
                             alpha_comparison_boxplot_fp,
                             _index_headers['alpha_diversity']))
                else:
                    logger.write("Skipping compare_alpha_diversity.py"
                                 " for %s as %s exists.\n\n"
                                 % (alpha_metric, compare_alpha_output_dir))
        else:
            logger.write("Skipping compare_alpha_diversity.py as"
                         " no categories were provided.\n\n")

    if not suppress_taxa_summary:
        taxa_plots_output_dir = '%s/taxa_plots/' % output_dir
        # need to check for existence of any html files, since the user can
        # select only certain ones to be generated
        existing_taxa_plot_html_fps = glob(join(taxa_plots_output_dir,
                                                'taxa_summary_plots', '*.html'))
        if len(existing_taxa_plot_html_fps) == 0:
            run_summarize_taxa_through_plots(
                otu_table_fp=biom_fp,
                mapping_fp=mapping_fp,
                output_dir=taxa_plots_output_dir,
                mapping_cat=None,
                sort=True,
                command_handler=command_handler,
                params=params,
                qiime_config=qiime_config,
                logger=logger,
                suppress_md5=True,
                status_update_callback=status_update_callback)
        else:
            logger.write("Skipping summarize_taxa_through_plots.py for as %s exist(s).\n\n"
                         % ', '.join(existing_taxa_plot_html_fps))

        index_links.append(('Taxa summary bar plots',
                            '%s/taxa_summary_plots/bar_charts.html'
                            % taxa_plots_output_dir,
                            _index_headers['taxa_summary']))
        index_links.append(('Taxa summary area plots',
                            '%s/taxa_summary_plots/area_charts.html'
                            % taxa_plots_output_dir,
                            _index_headers['taxa_summary']))
        for category in categories:
            taxa_plots_output_dir = '%s/taxa_plots_%s/' % (output_dir,
                                                           category)
            # need to check for existence of any html files, since the user can
            # select only certain ones to be generated
            existing_taxa_plot_html_fps = glob(
                '%s/taxa_summary_plots/*.html' %
                taxa_plots_output_dir)
            if len(existing_taxa_plot_html_fps) == 0:
                run_summarize_taxa_through_plots(
                    otu_table_fp=biom_fp,
                    mapping_fp=mapping_fp,
                    output_dir=taxa_plots_output_dir,
                    mapping_cat=category,
                    sort=True,
                    command_handler=command_handler,
                    params=params,
                    qiime_config=qiime_config,
                    logger=logger,
                    suppress_md5=True,
                    status_update_callback=status_update_callback)
            else:
                logger.write("Skipping summarize_taxa_through_plots.py for %s as %s exist(s).\n\n"
                             % (category, ', '.join(existing_taxa_plot_html_fps)))

            index_links.append(('Taxa summary bar plots',
                                '%s/taxa_summary_plots/bar_charts.html'
                                % taxa_plots_output_dir,
                                _index_headers['taxa_summary_categorical'] % category))
            index_links.append(('Taxa summary area plots',
                                '%s/taxa_summary_plots/area_charts.html'
                                % taxa_plots_output_dir,
                                _index_headers['taxa_summary_categorical'] % category))

    if not suppress_group_significance:
        params_str = get_params_str(params['group_significance'])
        # group significance tests, aka category significance
        for category in categories:
            group_signifance_fp = \
                '%s/group_significance_%s.txt' % (output_dir, category)
            if not exists(group_signifance_fp):
                # Build the OTU cateogry significance command
                group_significance_cmd = \
                    'group_significance.py -i %s -m %s -c %s -o %s %s' %\
                    (rarefied_biom_fp, mapping_fp, category,
                     group_signifance_fp, params_str)
                commands.append([('Group significance (%s)' % category,
                                  group_significance_cmd)])
            else:
                logger.write("Skipping group_significance.py for %s as %s exists.\n\n"
                             % (category, group_signifance_fp))

            index_links.append(('Category significance (%s)' % category,
                                group_signifance_fp,
                                _index_headers['group_significance']))

    filtered_biom_gzip_fp = '%s.gz' % filtered_biom_fp
    if not exists(filtered_biom_gzip_fp):
        commands.append(
            [('Compress the filtered BIOM table', 'gzip %s' %
              filtered_biom_fp)])
    else:
        logger.write("Skipping compressing of filtered BIOM table as %s exists.\n\n"
                     % filtered_biom_gzip_fp)
    index_links.append(
        ('Filtered BIOM table (minimum sequence count: %d)' % sampling_depth,
         filtered_biom_gzip_fp,
         _index_headers['run_summary']))

    rarefied_biom_gzip_fp = '%s.gz' % rarefied_biom_fp
    if not exists(rarefied_biom_gzip_fp):
        commands.append(
            [('Compress the rarefied BIOM table', 'gzip %s' %
              rarefied_biom_fp)])
    else:
        logger.write("Skipping compressing of rarefied BIOM table as %s exists.\n\n"
                     % rarefied_biom_gzip_fp)
    index_links.append(
        ('rarefied BIOM table (sampling depth: %d)' % sampling_depth,
         rarefied_biom_gzip_fp,
         _index_headers['run_summary']))

    if len(commands) > 0:
        command_handler(commands, status_update_callback, logger)
    else:
        logger.close()

    generate_index_page(index_links, index_fp)
