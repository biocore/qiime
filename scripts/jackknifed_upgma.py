#!/usr/bin/env python
# File created on 19 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.92"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from optparse import make_option
from os import makedirs
from qiime.util import load_qiime_config, parse_command_line_parameters,\
 raise_error_on_parallel_unavailable
from qiime.parse import parse_qiime_parameters
from qiime.workflow import run_jackknifed_upgma_clustering, print_commands,\
    call_commands_serially, print_to_stdout, no_status_updates

script_info={}
script_info['brief_description']="""A workflow script for performing jackknifed UPGMA clustering"""
script_info['script_description']="""To directly measure the robustness of individual UPGMA clusters, one can perform jackknifing (repeatedly resampling a subset of the available data from each sample). This process will utilize the following script: beta_diversity.py, multiple_rarefactions.py, upgma_cluster.py and tree_compare.py)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""These steps are performed by the following command:

1. Compute beta diversity distance matrix from otu table (and tree, if applicable)

2. Build rarefied OTU tables;

3. Build UPGMA tree from full distance matrix;

4. Compute distance matrics for rarefied OTU tables; 

5. Build UPGMA trees from rarefied OTU table distance matrices;

6. Compare rarefied OTU table distance matrix UPGMA trees to tree full UPGMA tree and write support file and newick tree with support values as node labels.

""","""jackknifed_upgma.py -i inseqs1_otu_table.txt -t inseqs1_rep_set.tre -p custom_parameters_jack.txt -o wf_jack -e 5 -v"""))
script_info['output_description']="""This scripts results in several distance matrices (from beta_diversity.py), several rarified otu tables (from multiple_rarefactions.py) several UPGMA trees (from upgma_cluster.py) and a supporting file and newick tree with support values (from tree_compare.py)."""

qiime_config = load_qiime_config()

script_info['required_options']=[\
 make_option('-i','--otu_table_fp',\
            help='the input fasta file [REQUIRED]'),\
 make_option('-o','--output_dir',\
            help='the output directory [REQUIRED]'),\
 make_option('-p','--parameter_fp',\
            help='path to the parameter file [REQUIRED]'),\
 make_option('-e','--seqs_per_sample',type='int',\
     help='number of sequences to include in each jackknifed subset'+\
            ' [REQUIRED]')
]

script_info['optional_options']=[\
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

script_info['version'] = __version__


def main():
    option_parser, opts, args = \
        parse_command_line_parameters(**script_info)
      
    verbose = opts.verbose
    
    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    tree_fp = opts.tree_fp
    seqs_per_sample = opts.seqs_per_sample
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
    
    run_jackknifed_upgma_clustering(otu_table_fp=otu_table_fp,\
     tree_fp=tree_fp,\
     seqs_per_sample=seqs_per_sample,\
     output_dir=output_dir, 
     command_handler=command_handler,
     params=parse_qiime_parameters(parameter_f),\
     qiime_config=qiime_config,\
     parallel=parallel,\
     status_update_callback=status_update_callback)


if __name__ == "__main__":
    main()
