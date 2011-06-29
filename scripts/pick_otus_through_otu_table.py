#!/usr/bin/env python
# File created on 27 Oct 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from qiime.util import make_option
from os import makedirs
from qiime.util import (load_qiime_config, 
                        parse_command_line_parameters,
                        get_options_lookup)
from qiime.parse import parse_qiime_parameters
from qiime.workflow import (run_qiime_data_preparation, print_commands,
    call_commands_serially, print_to_stdout, no_status_updates,
    validate_and_set_jobs_to_start)

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description'] = """A workflow script for picking OTUs through building OTU tables"""
script_info['script_description'] = """This script takes a sequence file and performs all processing steps through building the OTU table."""

script_info['script_usage'] = []

script_info['script_usage'].append(("""Simple example""","""The following command will start an analysis on inseq1.fasta (-i), which is a post-split_libraries fasta file. The sequence identifiers in this file should be of the form <sample_id>_<unique_seq_id>. The following steps, corresponding to the preliminary data preparation, are applied.

1. Pick OTUs with uclust at similarity of 0.97;

2. Pick a representative set with the most_abundant method;

3. Align the representative set with PyNAST;

4. Assign taxonomy with RDP classifier;

5. Filter the alignment prior to tree building - remove positions which are all gaps, and specified as 0 in the lanemask;

6. Build a phylogenetic tree with FastTree;

7. Build an OTU table.

All output files will be written to the directory specified by -o, and 
subdirectories as appropriate.
""","""pick_otus_through_otu_table.py -i inseqs1.fasta -o wf1/ -p custom_parameters.txt"""))

script_info['output_description'] ="""This script will produce an OTU mapping file (pick_otus.py), a representative set of sequences (FASTA file from pick_rep_set.py), a sequence alignment file (FASTA file from align_seqs.py), taxonomy assignment file (from assign_taxonomy.py), a filtered sequence alignment (from filter_alignment.py), a phylogenetic tree (Newick file from make_phylogeny.py) and an OTU table (from make_otu_table.py)."""

script_info['required_options'] = [
    make_option('-i','--input_fp',
        help='the input fasta file [REQUIRED]'),
    make_option('-o','--output_dir',
        help='the output directory [REQUIRED]'),
    ]

script_info['optional_options'] = [\
 make_option('-p','--parameter_fp',
    help='path to the parameter file, which specifies changes'+\
        ' to the default behavior. '+\
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .'+\
        ' [if omitted, default values will be used]'),
 make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
 make_option('-w','--print_only',action='store_true',\
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]',default=False),\
 make_option('-a','--parallel',action='store_true',\
        dest='parallel',default=False,\
        help='Run in parallel where available [default: %default]'),
 options_lookup['jobs_to_start_workflow']
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose
    
    input_fp = opts.input_fp
    output_dir = opts.output_dir
    verbose = opts.verbose
    print_only = opts.print_only
    
    parallel = opts.parallel
    # No longer checking that jobs_to_start > 2, but
    # commenting as we may change our minds about this.
    #if parallel: raise_error_on_parallel_unavailable()
    
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
    
    run_qiime_data_preparation(
     input_fp, 
     output_dir,
     command_handler=command_handler,
     params=params,
     qiime_config=qiime_config,
     parallel=parallel,\
     status_update_callback=status_update_callback)

if __name__ == "__main__":
    main()    
