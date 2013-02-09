#!/usr/bin/env python
# File created on 09 Feb 2013
from __future__ import division
import re
from os import makedirs, listdir
from glob import glob
from os.path import split, splitext, join, dirname, abspath
from datetime import datetime
from numpy import array
from cogent.util.misc import safe_md5
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import (parse_mapping_file, 
                        parse_qiime_parameters,
                        mapping_file_to_dict)
from qiime.util import (compute_seqs_per_library_stats,
                        get_qiime_scripts_dir,
                        create_dir, guess_even_sampling_depth,
                        get_interesting_mapping_fields,qiime_system_call,
                        get_qiime_library_version)
from biom.parse import parse_biom_table
from cogent.core.moltype import IUPAC_DNA_ambiguities
import os
from qiime.workflow import (print_to_stdout,
                            run_beta_diversity_through_plots,
                            run_qiime_alpha_rarefaction,
                            generate_log_fp,
                            WorkflowLogger,
                            log_input_md5s,
                            call_commands_serially,
                            get_params_str,
                            run_summarize_taxa_through_plots)

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

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


def run_core_qiime_analyses(
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
    status_update_callback=print_to_stdout):
    """
    """
    # Validate categories provided by the users
    mapping_data, mapping_categories, _ = parse_mapping_file(open(mapping_fp,'U'))
    if categories != None:
        for c in categories:
            try:
                cat_idx = mapping_categories.index(c)
            except IndexError:
                raise ValueError, ("Category '%s' is not a column header "
                 "in your mapping file. "
                 "Categories are case and white space sensitive. Valid "
                 "choices are: (%s)" % (c,', '.join(mapping_categories)))
            category_values = []
            for e in mapping_data:
                category_values.append(e[cat_idx])
            if len(set(category_values)) < 2:
                raise ValueError, ("Category '%s' contains only one value. "
                 "Categories analyzed here require at least two values.")
            
    else:
        categories= []
    
    # prep some variables
    if params == None:
        params = parse_qiime_parameters([])
        
    create_dir(output_dir)
    index_fp = '%s/index.html' % output_dir
    index_links = []
    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()
    
    # begin logging
    log_fp = generate_log_fp(output_dir)
    index_links.append(('Master run log',log_fp,'Log files'))
    logger = WorkflowLogger(log_fp,
                            params=params,
                            qiime_config=qiime_config)
    input_fps = [biom_fp,mapping_fp]
    if tree_fp != None:
        input_fps.append(tree_fp)
    log_input_md5s(logger,input_fps)
    
    
    bdiv_even_output_dir = '%s/bdiv_even%d/' % (output_dir,sampling_depth)
    even_dm_fps = run_beta_diversity_through_plots(
     otu_table_fp=biom_fp, 
     mapping_fp=mapping_fp,
     output_dir=bdiv_even_output_dir,
     command_handler=command_handler,
     params=params,
     qiime_config=qiime_config,
     sampling_depth=sampling_depth,
     histogram_categories=categories,
     tree_fp=tree_fp,
     parallel=parallel,
     logger=logger,
     status_update_callback=status_update_callback)
     
    for bdiv_metric, dm_fp in even_dm_fps:
        for category in categories:
            ## Add category-specific comparisons
            pass
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
    arare_full_output_dir = '%s/arare_max%d/' % (output_dir,sampling_depth)
    run_qiime_alpha_rarefaction(
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
     status_update_callback=status_update_callback)
    
    index_links.append(('Alpha rarefaction plots',
                        '%s/alpha_rarefaction_plots/rarefaction_plots.html'\
                          % arare_full_output_dir,
                        "Alpha rarefaction results"))
    
    taxa_plots_output_dir = '%s/taxa_plots/' % output_dir
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
     status_update_callback=status_update_callback)

    index_links.append(('Taxa summary bar plots',
                        '%s/taxa_summary_plots/bar_charts.html'\
                          % taxa_plots_output_dir,
                        "Taxonomic summary results"))
    index_links.append(('Taxa summary area plots',
                        '%s/taxa_summary_plots/area_charts.html'\
                          % taxa_plots_output_dir,
                        "Taxonomic summary results"))
    for c in categories:
        taxa_plots_output_dir = '%s/taxa_plots_%s/' % (output_dir,c)
        run_summarize_taxa_through_plots(
         otu_table_fp=biom_fp,
         mapping_fp=mapping_fp,
         output_dir=taxa_plots_output_dir,
         mapping_cat=c, 
         sort=True,
         command_handler=command_handler,
         params=params,
         qiime_config=qiime_config,
         logger=logger, 
         status_update_callback=status_update_callback)

        index_links.append(('Taxa summary bar plots',
                            '%s/taxa_summary_plots/bar_charts.html'\
                              % taxa_plots_output_dir,
                            "Taxonomic summary results (by %s)" % c))
        index_links.append(('Taxa summary area plots',
                            '%s/taxa_summary_plots/area_charts.html'\
                              % taxa_plots_output_dir,
                            "Taxonomic summary results (by %s)" % c))
    
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
         (biom_fp, mapping_fp, category, 
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
         (biom_fp, mapping_fp, category, 
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
