#!/usr/bin/env python
# File created on 29 May 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from os import makedirs
from qiime.util import load_qiime_config
from qiime.parse import parse_qiime_parameters
from qiime.workflow import (run_summarize_taxa_through_plots, print_commands,
    call_commands_serially, print_to_stdout, no_status_updates,
    validate_and_set_jobs_to_start)

qiime_config = load_qiime_config()

#summarize_taxa_through_plots.py
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""A workflow script for performing taxonomy summaries and plots"""
script_info['script_description']="""
The steps performed by this script are:

1. Summarize OTU by Category

2. Summarize Taxonomy

3. Plot Taxonomy Summary

"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples""","""""","""summarize_taxa_through_plots.py -o wf_taxa_sum -i otu_table.txt -m inseqs1_mapping.txt -p custom_parameters.txt"""))
script_info['script_usage'].append(("""""","""Alternatively, the user can supply a mapping_category, where the OTU is summarized based on a mapping category:""","""summarize_taxa_through_plots.py -o wf_taxa_sum -i otu_table.txt -m inseqs1_mapping.txt -p custom_parameters.txt -c Treatment"""))
script_info['output_description']="""The results of this script is a folder ("wf_taxa_sum/") containing taxonomy summary files (at different levels) and a folder containing taxonomy summary plots. Additionally, if a mapping_catgory is supplied there will be a summarized OTU table."""
script_info['required_options']=[\
    make_option('-i','--otu_table_fp',\
        help='the input otu table [REQUIRED]'),\
    make_option('-m','--mapping_fp',\
        help='path to the mapping file [REQUIRED]'),\
    make_option('-o','--output_dir',\
        help='the output directory [REQUIRED]'),
]
script_info['optional_options']=[\
    make_option('-p','--parameter_fp',
        help='path to the parameter file, which specifies changes'+\
        ' to the default behavior. '+\
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters.'+\
        ' [if omitted, default values will be used]'),
    make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
    make_option('-w','--print_only',action='store_true',\
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]',default=False),\
    make_option('-c','--mapping_category',default=None,
        help='Summarize OTU table using this category. [default: %default]'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    verbose = opts.verbose
    
    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    mapping_fp = opts.mapping_fp
    verbose = opts.verbose
    print_only = opts.print_only
    mapping_cat = opts.mapping_category
    
    if opts.parameter_fp:
        try:
            parameter_f = open(opts.parameter_fp)
        except IOError:
            raise IOError,\
                "Can't open parameters file (%s). Does it exist? Do you have read access?" \
                 % opts.parameter_fp
        params = parse_qiime_parameters(parameter_f)
    else:
        params = parse_qiime_parameters([]) 
        # empty list returns empty defaultdict for now
                                   
    try:
        makedirs(output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            print "Output directory already exists. Please choose "+\
             "a different directory, or force overwrite with -f."
            exit(1)
        
    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
    
    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates
     
    run_summarize_taxa_through_plots(
     otu_table_fp=otu_table_fp,
     mapping_fp=mapping_fp,
     output_dir=output_dir,
     mapping_cat=mapping_cat,
     command_handler=command_handler,
     params=params,
     qiime_config=qiime_config,
     status_update_callback=status_update_callback)

if __name__ == "__main__":
    main()