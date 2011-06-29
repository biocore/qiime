#!/usr/bin/env python
# File created on 04 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from os import makedirs
from qiime.util import load_qiime_config
from qiime.parse import parse_qiime_parameters
from qiime.workflow import (run_qiime_alpha_rarefaction, print_commands,
    call_commands_serially, print_to_stdout, no_status_updates,
    validate_and_set_jobs_to_start)

qiime_config = load_qiime_config()

#alpha_rarefaction.py
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""A workflow script for performing alpha rarefaction"""
script_info['script_description']="""
The steps performed by this script are:

1. Generate rarefied OTU tables;

2. Compute alpha diversity metrics for each rarefied OTU table;

3. Collate alpha diversity results;

4. Generate alpha rarefaction plots.
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example""","""""","""alpha_rarefaction.py -o rare1 -i otu_table.txt -t inseqs1_rep_set.tre -m inseqs1_mapping.txt -p custom_parameters.txt"""))
script_info['output_description']="""The results of this script is a folder ("rare1/") containing rarefied otu tables, alpha-diversity for each otu table, a file containing the results from collating the alpha-diversity results and a folder containing the rarefaction plots."""
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
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .'+\
        ' [if omitted, default values will be used]'),
 make_option('-n','--num_steps',type='int',\
     help='number of steps (or rarefied OTU table sizes) to make between '+\
      'min and max counts [default: %default]',default=10),\
 make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
 make_option('-w','--print_only',action='store_true',\
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]',default=False),\
 make_option('-a','--parallel',action='store_true',\
        dest='parallel',default=False,\
        help='Run in parallel where available [default: %default]'),\
 make_option('-t','--tree_fp',\
            help='path to the tree file [default: %default; '+\
            'REQUIRED for phylogenetic measures]'),
 options_lookup['jobs_to_start_workflow']]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    verbose = opts.verbose
    
    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    mapping_fp = opts.mapping_fp
    tree_fp = opts.tree_fp
    num_steps = opts.num_steps
    verbose = opts.verbose
    print_only = opts.print_only
    parallel = opts.parallel
    
    if opts.parameter_fp:
        try:
            parameter_f = open(opts.parameter_fp)
        except IOError:
            raise IOError,\
             "Can't open parameters file (%s). Does it exist? Do you have read access?"\
             % opts.parameter_fp
        params = parse_qiime_parameters(parameter_f)
    else:
        params = parse_qiime_parameters([]) 
        # empty list returns empty defaultdict for now
    
    jobs_to_start = opts.jobs_to_start
    default_jobs_to_start = qiime_config['jobs_to_start']
    validate_and_set_jobs_to_start(params,
                                   jobs_to_start,
                                   default_jobs_to_start,
                                   parallel,
                                   option_parser)
                                   
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
     
    run_qiime_alpha_rarefaction(otu_table_fp=otu_table_fp,\
     mapping_fp=mapping_fp,\
     output_dir=output_dir,\
     command_handler=command_handler,\
     params=params,
     qiime_config=qiime_config,\
     tree_fp=tree_fp,\
     num_steps=num_steps,\
     parallel=parallel,\
     status_update_callback=status_update_callback)

if __name__ == "__main__":
    main()
