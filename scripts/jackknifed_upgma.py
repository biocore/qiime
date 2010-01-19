#!/usr/bin/env python
# File created on 19 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, Qiime"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"


from optparse import OptionParser
from os import makedirs
from qiime.util import load_qiime_config
from qiime.parse import parse_qiime_parameters
from qiime.workflow import run_jackknifed_upgma_clustering, print_commands,\
    call_commands_serially, print_to_stdout, no_status_updates

usage_str = """usage: %prog [options] {-i OTU_TABLE_FP -o OUTPUT_DIR -p PARAMETERS_FP -e SEQS_PER_SAMPLE}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

 The steps performed by this script are:
  1) Compute beta diversity distance matrix from otu table (and
   tree, if applicable)
  2) Build rarefied OTU tables;
  3) Build UPGMA tree from full distance matrix;
  4) Compute distance matrics for rarefied OTU tables;
  5) Build UPGMA trees from rarefied OTU table distance matrices;
  6) Compare rarefied OTU table distance matrix UPGMA trees 
   to tree full UPGMA tree and write support file and newick tree
   with support values as node labels.

Example usage:
 python ~/code/Qiime/scripts/jackknifed_upgma.py -i inseqs1_otu_table.txt -t inseqs1_rep_set.tre -p custom_parameters_jack.txt -o wf_jack -e 5 -v
"""

qiime_config = load_qiime_config()

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--otu_table_fp',\
            help='the input fasta file [REQUIRED]')
    parser.add_option('-o','--output_dir',\
            help='the output directory [REQUIRED]')
    parser.add_option('-p','--parameter_fp',\
            help='path to the parameter file [REQUIRED]')
    parser.add_option('-e','--seqs_per_sample',type='int',\
     help='number of sequences to include in each jackknifed subset'+\
            ' [REQUIRED]')
    parser.add_option('-t','--tree_fp',\
            help='path to the tree file [default: %default; '+\
            'REQUIRED for phylogenetic measures]')

    parser.add_option('-v','--verbose',action='store_true',\
        dest='verbose',help='Print information during execution -- '+\
        'useful for debugging [default: %default]')
        
    parser.add_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]')
        
    parser.add_option('-w','--print_only',action='store_true',\
        dest='print_only',help='Print the commands but don\'t call them -- '+\
        'useful for debugging [default: %default]')
        
    parser.add_option('-s','--serial',action='store_true',\
        dest='serial',help='Do not use parallel scripts'+\
        ' [default: %default]')

    parser.set_defaults(verbose=False,print_only=False,\
     serial=False)

    required_options = ['otu_table_fp', 'output_dir','parameter_fp',\
     'seqs_per_sample']

    opts,args = parser.parse_args()
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option)

    return opts,args


if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    
    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    tree_fp = opts.tree_fp
    seqs_per_sample = opts.seqs_per_sample
    verbose = opts.verbose
    print_only = opts.print_only
    
    parallel = not opts.serial
    if parallel:
        # Keeping this check in until fully tested
        print "Parallel runs not yet supported. Running in single proc mode."
        parallel = False
    
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