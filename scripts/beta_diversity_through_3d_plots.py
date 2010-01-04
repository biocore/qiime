#!/usr/bin/env python
# File created on 04 Jan 2010.
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
from qiime.workflow import run_beta_diversity_through_3d_plot, print_commands,\
    call_commands_serially, print_to_stdout, no_status_updates
    
usage_str = """usage: %prog [options] {-i OTU_TABLE_FP -m MAPPING_FP -o OUTPUT_DIR}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

REQUIRED:
 You must edit the following parameters in a custom parameters file:
  beta_diversity:metric
  
 This is the value that would be passed to beta_diversity.py via -m/--metric.
 
 The steps performed by this script are:
  1) Compute a beta diversity distance matrix;
  2) Peform a principle coordinates analysis on the result of Step 1;
  3) Generate a 3D prefs file for optimized coloring of continuous variables;
  4) Generate a 3D plot for all mapping fields with colors optimized for 
   continuous data;
  5) Generate a 3D plot for all mapping fields with colors optimized for 
   discrete data.
   
  This script is in early development status.

Example usage:

 python beta_diversity_through_3d_plots.py -i otu_table.txt -o bdiv1 -t inseqs1_rep_set.tre -m inseqs1_mapping.txt -p custom_parameters.txt
"""

qiime_config = load_qiime_config()

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--otu_table_fp',\
            help='the input fasta file [REQUIRED]')
    parser.add_option('-m','--mapping_fp',\
            help='path to the mapping file [REQUIRED]')
    parser.add_option('-o','--output_dir',\
            help='the output directory [REQUIRED]')
    parser.add_option('-t','--tree_fp',\
            help='path to the tree file [default: %default; '+\
            'REQUIRED for phylogenetic measures]')
            
    parser.add_option('-p','--parameter_fp',\
     help='path to the parameter file '+\
          '[default: %default]')

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
    
    # parallel not supported yet!    
    #   SHOULD MAKE THIS -s FOR SERIAL/SINGLE PROC MODE -- PARALLEL WILL
    #   BE THE DEFAULT ONCE IT'S READY TO GO.
    # parser.add_option('-p','--parallel',action='store_true',\
    #     dest='parallel',help='Use parallel scripts where applicable'+\
    #     ' [default: %default]')

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False,print_only=False,parallel=False,
     parameter_fp=qiime_config['qiime_parameter_fp'])

    required_options = ['otu_table_fp', 'output_dir', 'mapping_fp']

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
    mapping_fp = opts.mapping_fp
    tree_fp = opts.tree_fp
    verbose = opts.verbose
    print_only = opts.print_only
    
    ## REMEMBER TO GRAB opts.parallel WHEN IT BECOMES SUPPORTED.
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
     
    run_beta_diversity_through_3d_plot(otu_table_fp=otu_table_fp,\
     mapping_fp=mapping_fp,\
     output_dir=output_dir,\
     command_handler=command_handler,\
     params=parse_qiime_parameters(parameter_f),\
     qiime_config=qiime_config,\
     tree_fp=tree_fp,\
     parallel=parallel,\
     status_update_callback=status_update_callback)