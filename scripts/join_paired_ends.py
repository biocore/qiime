#!/usr/bin/env python
# file: join_paired_ends.py

__author__ = "Mike Robeson"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Mike Robeson"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Mike Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"

from cogent.parse.fastq import MinimalFastqParser
from qiime.join_paired_ends import (join_method_names,
                                    join_method_constructors,
                                    write_synced_barcodes_fastq)
from qiime.util import (parse_command_line_parameters, get_options_lookup, 
                        make_option, load_qiime_config, create_dir)
import os
import gzip


options_lookup = get_options_lookup()
qiime_config = load_qiime_config()

script_info={}
script_info['brief_description']= """Joins paired-end Illumina reads."""
script_info['script_description'] = """This script take forward and reverse Illumina reads and joins them based on the method chosen. Will optionally, create an updated index reads file to match the surviving joined paired-ends. If the option to write an updated index file is chosen, be sure that the order of the index-reads file is in the same order as the reads files you plan to join.

Currently, there are two methods that can be selected by the user to join paired-end data:

1. fastq-join - Erik Aronesty, 2011. ea-utils : "Command-line tools for processing biological sequencing data" (http://code.google.com/p/ea-utils)

2. SeqPrep - (https://github.com/jstjohn/SeqPrep)
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Join paired-ends with \'fastq-join\':""","""This is the default method to join paired-end Illumina data:""",""" %prog -f $PWD/forward_reads.fastq -r reverse_reads.fastq"""))
script_info['script_usage'].append(("""Join paired-ends with \'SeqPrep\':""","""Produces similar output to the \'fastq-join\' but returns data in gzipped format.""",""" %prog -m SeqPrep -f $PWD/forward_reads.fastq -r reverse_reads.fastq"""))
script_info['script_usage'].append(("""Update the index / barcode reads file to match the surviving joined pairs.""","""This is required if you will be using \'split_libraries_fastq.py\'.""",""" %prog -f $PWD/forward_reads.fastq -r reverse_reads.fastq -b index.reads.fastq"""))
script_info['output_description'] = """All paired-end joining software will return a joined / merged / assembled paired-end fastq file. Depending on the method chosen, additional files may be written to the user-specified output directory. 


1. fastq-join will output fastq-formatted files as:
   \"*.join\" - assembled / joined reads output
   \"*.un1\" - unassembled / unjoined reads1 output
   \"*.un2\" - unassembled / unjoined reads2 output

2. SeqPrep will output fastq-formatted gzipped files as: 
   \"*_assembled.gz\" - unassembled / unjoined reads1 output
   \"*_unassembled_R1.gz\" - unassembled / unjoined reads1 output
   \"*_unassembled_R2.gz\" - unassembled / unjoined reads2 output

3. If a barcode / index file is provided via the \'-b\' option, an updated
   barcodes file will be output as:
   \"..._barcodes.fastq\"
    This barcode / index file must be used in conjunction with the joined
    paired-ends file as input to \'split_libraries_fastq.py\'. Except for
    missing reads that may result from failed merging of paired-ends, the
    index-reads and joined-reads must be in the same order.

"""
script_info['required_options'] = [\
    make_option('-f','--forward_reads_fp',type="existing_filepath",
                 dest='forward_reads_fp',
                 help='Path to read input forward reads in FASTQ format.'),
    make_option('-r','--reverse_reads_fp',type="existing_filepath",
                 dest='reverse_reads_fp',
                 help='Path to read input reverse reads in FASTQ format.')]
script_info['optional_options'] = [\
    make_option('-m', '--pe_join_method', action='store', type='choice',
                choices=list(join_method_names.keys()),
                help='Method to use for joining paired-ends. Valid choices'+\
                      ' are: ' + ', '.join(join_method_names.keys())+\
                      ' [default: %default]', default='fastq-join'),
    make_option('-o', '--output_dir', action='store', type='new_dirpath',\
                help='Path to store '+\
                      'result file [default: <PE_JOIN_METHOD>_joined]'),
    make_option('-b','--index_reads_fp',type='existing_filepath',
                dest='index_reads_fp',
                help='Path to read the barcode / index reads in FASTQ format.'
                ' Will be filtered based on surviving joined pairs.'),
    make_option('-j', '--min_overlap', 
                help='Applies to both fastq-join and SeqPrep methods.'+\
                      ' Minimum allowed overlap in base-pairs required to join pairs.'+\
                      ' Defaults to recomended settings: fastq-join (6), SeqPrep (15) '+\
                      ', [default: %default]', default='default'),
    make_option('-p', '--perc_max_diff', action='store', type='int',
                help='Only applies to fastq-join method, otherwise ignored.'+\
                     ' Maximum allowed % differences within region of overlap'+\
                      ',  [default: %default]', default=8),
    make_option('-y', '--max_ascii_score', action='store', type='string',
                help='Only applies to SeqPrep method, otherwise ignored.'+\
                      ' Maximum quality score / ascii code allowed to appear within'+\
                      ' joined pairs output. For more information see:'+\
                      ' http://en.wikipedia.org/wiki/FASTQ_format ' 
                      ' [default: %default]', default='J'),
    make_option('-n', '--min_frac_match', action='store', type='float',
                help='Only applies to SeqPrep method, otherwise ignored.'+\
                      ' Minimum allowed fraction of matching bases required to join reads'+\
                      ',  [default: %default]', default=0.9),
    make_option('-g', '--max_good_mismatch', action='store', type='float',
                help='Only applies to SeqPrep method, otherwise ignored.'+\
                      ' Maximum mis-matched high quality bases allowed'+\
                      ' to join reads. ' + '[default: %default]', default=0.2),
    make_option('-6', '--phred_64', action='store_true',
                help='Only applies to SeqPrep method, otherwise ignored.'+\
                      ' Set if input reads are in phred+64 format. Output will '\
                      'always be phred+33. [default: %default]',
                        default=False)]

script_info['version'] = __version__


def main():
    # parse command line parameters
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    # Create local copy of options
    forward_reads_fp = opts.forward_reads_fp 
    reverse_reads_fp = opts.reverse_reads_fp
    pe_join_method = opts.pe_join_method
    output_dir = opts.output_dir
    # fastq-join only options:
    perc_max_diff = opts.perc_max_diff
    # SeqPrep only options:
    max_ascii_score = opts.max_ascii_score
    min_frac_match = opts.min_frac_match
    max_good_mismatch = opts.max_good_mismatch
    phred_64 = opts.phred_64
    # both fastq-join & SeqPrep options
    min_overlap = opts.min_overlap


    # check for valid paired-end join method:
    if not (pe_join_method in join_method_constructors or 
            pe_join_method in join_method_names):
       option_parser.error(\
        'Invalid paired-end join method: %s. \nValid choces are: %s'\
        %opts.pe_join_method,' '.join(join_method_constructors.keys() +
                             join_method_names.keys()))

    # set default min_overlap values according to join method
    if min_overlap != "default":
        try:
            min_overlap = int(min_overlap)
        except ValueError:
            raise ValueError, ("--min_overlap must either be 'default'", 
               "or an int value")
    if min_overlap == "default":
        if pe_join_method == "fastq-join":
            min_overlap = 6
        elif pe_join_method == "SeqPrep":
            min_overlap = 15

    # determine output directory:   
    if output_dir: # user specified output directory
        output_dir = os.path.abspath(opts.output_dir)
    else: # default output dir to location of infile
        output_dir = os.path.join(os.path.dirname(os.path.abspath(
                                  forward_reads_fp)),
                                  pe_join_method + '_joined')
    
    create_dir(output_dir, fail_on_exist=False)

    # send parameters to appropriate join method
    # currently only two join methods exist:
    # 'fastq-join' and 'SeqPrep'
    if pe_join_method == "fastq-join":
        join_func = join_method_names["fastq-join"]
        paths = join_func(forward_reads_fp,
                          reverse_reads_fp,
                          perc_max_diff = perc_max_diff,
                          min_overlap = min_overlap,
                          working_dir = output_dir)

    if pe_join_method == "SeqPrep":
        join_func = join_method_names["SeqPrep"]
        paths = join_func(forward_reads_fp,
                          reverse_reads_fp,
                          max_overlap_ascii_q_score = max_ascii_score,
                          min_overlap = min_overlap,
                          max_mismatch_good_frac = max_good_mismatch,
                          min_frac_matching = min_frac_match,
                          phred_64 = phred_64,
                          working_dir = output_dir)

   
    # If index / barcode file is supplied, filter unused barcode reads
    # and write them to a new file. Name based on joined-pairs / assembled
    # outfile 
    if opts.index_reads_fp:
        index_reads = opts.index_reads_fp
        assembly_fp = paths['Assembled'] # grab joined-pairs output path
        write_synced_barcodes_fastq(assembly_fp,index_reads)


if __name__ == "__main__":
    main()

