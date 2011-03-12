#!/usr/bin/env python
# File created on 27 Oct 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from optparse import make_option
from os import makedirs
from qiime.util import load_qiime_config, parse_command_line_parameters
from qiime.parse import parse_qiime_parameters
from qiime.workflow import run_core_qiime_analyses, print_commands,\
    call_commands_serially, print_to_stdout, no_status_updates

qiime_config = load_qiime_config()

script_info={}
script_info['brief_description'] = """A workflow script for running a core QIIME workflow."""
script_info['script_description'] = """This script plugs several QIIME steps together to form a basic full data analysis workflow. The steps include quality filtering and demultiplexing sequences, running the pick_otus_through_otu_table.py workflow (pick otus and representative sequences, assign taxonomy, align representative sequences, build a tree, and build and OTU table), generating 3d beta diversity PCoA plots, generating alpha rarefaction plots, identifying OTUs that are differentially represented in different categories, and several additional analysis."""
script_info['script_usage'] = [("","Run serial analysis","%prog -i Fasting_Example.fna -q Fasting_Example.qual -o FastingStudy -p custom_parameters.txt -m Fasting_Map.txt -c Treatment,DOB -e 100")]


script_info['output_description'] ="""
"""

script_info['required_options'] = [
    make_option('-i','--input_fnas',
        help='the input fna files (pre-split_libraries) [REQUIRED]'),
    make_option('-q','--input_quals',
        help='the input qual files (pre-split_libraries) [REQUIRED]'),
    make_option('-o','--output_dir',
        help='the output directory [REQUIRED]'),
    make_option('-p','--parameter_fp',
        help='path to the parameter file [REQUIRED]'),
    make_option('-m','--mapping_fp',
        help='the mapping filepath [REQUIRED]'),
    ]

script_info['optional_options'] = [\
 make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),
 make_option('-w','--print_only',action='store_true',
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]',default=False),
 make_option('-a','--parallel',action='store_true',
        dest='parallel',default=False,
        help='Run in parallel where available [default: %default]'),
 make_option('-s','--sff_fp',
        help='the sff file [REQUIRED for denoising]'),
 make_option('-e','--seqs_per_sample',type='int',
     help='depth of coverage for even sampling [default: %default]'),
 make_option('-t','--reference_tree_fp',
            help='path to the tree file if one should be used (otherwise de novo '+\
            ' tree will be used) [default: %default]'),
 make_option('-c','--categories',
            help='the categories to compare (for otu_category_significance,'+\
            'supervised_learning.py, and cluster_quality.py steps) '+\
            '[default: %default; skip these steps]'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose
    
    input_fnas = opts.input_fnas
    input_quals = opts.input_quals
    output_dir = opts.output_dir
    sampling_depth = opts.seqs_per_sample
    categories = opts.categories
    reference_tree_fp = opts.reference_tree_fp
    mapping_fp = opts.mapping_fp
    sff_fp = opts.sff_fp
    verbose = opts.verbose
    print_only = opts.print_only
    
    parallel = opts.parallel
    
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
    
    run_core_qiime_analyses(
        fna_fps=input_fnas,
        qual_fps=input_quals,
        mapping_fp=mapping_fp,
        output_dir=output_dir,
        command_handler=command_handler,
        params=parse_qiime_parameters(parameter_f),
        qiime_config=qiime_config,
        categories=categories,
        sampling_depth=sampling_depth,
        arare_min_seqs_per_sample=10,
        arare_num_steps=10,
        reference_tree_fp=reference_tree_fp,
        sff_input_fp=sff_fp,
        parallel=parallel,
        status_update_callback=status_update_callback)

if __name__ == "__main__":
    main()    
