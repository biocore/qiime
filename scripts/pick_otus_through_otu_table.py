#!/usr/bin/env python
# File created on 27 Oct 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

from optparse import OptionParser
from os import makedirs
from qiime.util import load_qiime_config
from qiime.parse import parse_qiime_parameters
from qiime.workflow import run_qiime_data_preparation, print_commands,\
    call_commands_serially, print_to_stdout, no_status_updates
    
usage_str = """usage: %prog [options] {-i INPUT_FP -o OUTPUT_DIR -p PARAMETER_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

WARNING: THE INTERFACE FOR THE WORKFLOW SCRIPT WAS CHANGED SIGNIFICANTLY
ON 29 DEC 2009. REVIEW THE NEW Qiime/qiime_parameters.txt FILE TO SEE HOW
TO SPECIFY PARAMETERS FOR THE INDIVIDUAL STEPS. COPY THAT FILE AND EDIT THE
COPY TO DEFINE CUSTOM ANALYSES.

REQUIRED:
 You must add values for the following parameters in a custom parameters file:
 align_seqs:template_fp
 filter_alignment:lane_mask_fp 
 
 These are the values that you would typically pass as --template_fp to 
  align_seqs.py and lane_mask_fp to filter_alignment.py, respectively.

Example usage:

The following command will start an analysis on inseq1.fasta (-i), which is a 
 post-split_libraries fasta file. The sequence identifiers in this file
 should be of the form <sample_id>_<unique_seq_id>. The following steps,
 corresponding to the preliminary data preparation, are applied:
  1) Pick OTUs with cdhit at similarity of 0.97;
  2) Pick a representative set with the most_abundant method;
  3) Align the representative set with PyNAST (REQUIRED: SPECIFY TEMPLATE 
   ALIGNMENT with align_seqs:template_fp in the parameters file);
  4) Assign taxonomy with RDP classifier;
  5) Filter the alignment prior to tree building - remove positions which 
   are all gaps, and specified as 0 in the lanemask (REQUIRED: SPECIFY LANEMASK 
   with filter_alignment:lane_mask_fp in the parameters file);
  6) Build a phylogenetic tree with FastTree;
  7) Build an OTU table.
 
 All output files will be written to the directory specified by -o, and 
 subdirectories as appropriate.

pick_otus_through_otu_table.py -i /Users/caporaso/data/qiime_test_data/workflow/inseqs1.fasta -o /Users/caporaso/data/qiime_test_data/workflow/wf1/ -p /Users/caporaso/data/qiime_test_data/workflow/custom_parameters.txt
"""
qiime_config = load_qiime_config()

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--input_fp',\
            help='the input fasta file [REQUIRED]')
    parser.add_option('-o','--output_dir',\
            help='the output directory [REQUIRED]')
    parser.add_option('-p','--parameter_fp',\
            help='path to the parameter file [REQUIRED]')

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

    parser.set_defaults(verbose=False,print_only=False,serial=False)

    required_options = ['input_fp', 'output_dir', 'parameter_fp']

    opts,args = parser.parse_args()
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option)

    return opts,args


if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    
    input_fp = opts.input_fp
    output_dir = opts.output_dir
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
    
    run_qiime_data_preparation(input_fp, output_dir,\
     command_handler=command_handler,\
     params=parse_qiime_parameters(parameter_f),\
     qiime_config=qiime_config,\
     parallel=parallel,\
     status_update_callback=status_update_callback)
    