#!/usr/bin/env python
# File created on 04 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

from optparse import make_option
from os import makedirs
from qiime.util import load_qiime_config, parse_command_line_parameters,\
 raise_error_on_parallel_unavailable
from qiime.parse import parse_qiime_parameters
from qiime.workflow import run_beta_diversity_through_3d_plot, print_commands,\
    call_commands_serially, print_to_stdout, no_status_updates

script_description = """A workflow script for computing beta diversity distance matrices and the corresponding 3D plots

REQUIRED:
 You must edit the following parameters in a custom parameters file:
  beta_diversity:metric
  
 This is the value that would be passed to beta_diversity.py via -m/--metric."""

script_usage = """ The following steps are performed by the command below:
  1) Compute a beta diversity distance matrix;
  2) Peform a principle coordinates analysis on the result of Step 1;
  3) Generate a 3D prefs file for optimized coloring of continuous variables;
  4) Generate a 3D plot for all mapping fields with colors optimized for 
   continuous data;
  5) Generate a 3D plot for all mapping fields with colors optimized for 
   discrete data.
python beta_diversity_through_3d_plots.py -i otu_table.txt -o bdiv1 -t inseqs1_rep_set.tre -m inseqs1_mapping.txt -p custom_parameters.txt"""

qiime_config = load_qiime_config()

required_options = [\
 make_option('-i','--otu_table_fp',\
            help='the input fasta file [REQUIRED]'),\
 make_option('-m','--mapping_fp',\
            help='path to the mapping file [REQUIRED]'),\
 make_option('-o','--output_dir',\
            help='the output directory [REQUIRED]'),\
 make_option('-p','--parameter_fp',\
            help='path to the parameter file [REQUIRED]')
]

optional_options = [\
 make_option('-t','--tree_fp',\
            help='path to the tree file [default: %default; '+\
            'REQUIRED for phylogenetic measures]'),\
 make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
 make_option('-w','--print_only',action='store_true',\
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]',default=False),\
 make_option('-a','--parallel',action='store_true',\
        dest='parallel',default=False,\
        help='Run in parallel where available [default: %default]')
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    verbose = opts.verbose
    
    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    mapping_fp = opts.mapping_fp
    tree_fp = opts.tree_fp
    verbose = opts.verbose
    print_only = opts.print_only
    
    parallel = opts.parallel
    if parallel: raise_error_on_parallel_unavailable()
    
    try:
        parameter_f = open(opts.parameter_fp)
    except IOError:
        raise IOError,\
         "Can't open parameters file (%s). Does it exist? Do you have read access?"\
         % opts.parameter_fp
    
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
     
    run_beta_diversity_through_3d_plot(otu_table_fp=otu_table_fp,\
     mapping_fp=mapping_fp,\
     output_dir=output_dir,\
     command_handler=command_handler,\
     params=parse_qiime_parameters(parameter_f),\
     qiime_config=qiime_config,\
     tree_fp=tree_fp,\
     parallel=parallel,\
     status_update_callback=status_update_callback)

if __name__ == "__main__":
    main()