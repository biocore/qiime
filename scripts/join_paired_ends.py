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
from os.path import abspath, dirname, splitext
import gzip

# TODO:
# add method specific options
# rename output files to something more consistent between methods

options_lookup = get_options_lookup()
qiime_config = load_qiime_config()

script_info={}
script_info['brief_description']= """Joins paired-end Illumina reads."""
script_info['script_description'] = """This script take forward and reverse Illumina reads and joins them based on the method chosen. Will optionally, create an updated index reads file to match the surviving joined paired-ends. Currently, there are two methods that can be selected by the user:

1. fastq-join - Erik Aronesty, 2011. ea-utils : "Command-line tools for processing biological sequencing data" (http://code.google.com/p/ea-utils)

2. SeqPrep - (https://github.com/jstjohn/SeqPrep)
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Join paired-ends with \'fastq-join\':""","""This is the default method to join paired-end Illumina data:""",""" %prog -f $PWD/forward_reads.fastq -r reverse_reads.fastq"""))
script_info['script_usage'].append(("""Join paired-ends with \'SeqPrep\':""","""Produces similar output to the \'fastq-join\' but returns data in gzipped format.""",""" %prog -m SeqPrep -f $PWD/forward_reads.fastq -r reverse_reads.fastq"""))
script_info['script_usage'].append(("""Update the index / barcode reads file to match the surviving joined pairs.""","""This is required if you will be using \'split_libraries_fastq.py\'.""",""" %prog -f $PWD/forward_reads.fastq -r reverse_reads.fastq -b index.reads.fastq"""))
script_info['output_description'] = """All paired-end joining software will return a joined / merged / assembled paired-end fastq file. Depending on the method chosen, additional files may be written to the user-specified output directory. 


1. fastq-join will output fastq-formatted files as:
   \"..._join\" - assembled / joined reads output
   \"..._un1\" - unassembled / unjoined reads1 output
   \"..._un2\" - unassembled / unjoined reads2 output

2. SeqPrep will output fastq-formatted gzipped files as: 
   \"..._assembled.gz\" - unassembled / unjoined reads1 output
   \"..._unassembled_R1.gz\" - unassembled / unjoined reads1 output
   \"..._unassembled_R2.gz\" - unassembled / unjoined reads2 output

3. If a barcode / index file is provided via the \'-b\' option, an updated
   barcodes file will be output as:
   \"..._barcodes.fastq\"
    This barcode / index file must be used in conjunction with the joined
    paired-ends file as input to \'split_libraries_fastq.py\'.

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
                      'result file [default: <JOINED_METHOD>_joined]'),
    make_option('-b','--index_reads_fp',type='existing_filepath',
                dest='index_reads_fp',
                help='Path to read the barcode / index reads in FASTQ format.'
                ' Will be filtered based on surviving joined pairs.'),
    make_option('-p', '--perc_max_diff', action='store', type='int',
                help='fastq-join option:'+\
                     ' Maximum allowed % differences within region of overlap'+\
                      ',  [default: %default]', default=8),
    make_option('-j', '--min_overlap', action='store', type='int',
                help='fastq-join and SeqPrep option:'+\
                      ' Minimum allowed overlap in base-pairs required to join pairs.'+\
                      ' Recomended default settings are: fastq-join (6), SeqPrep (15) '+\
                      ', [default: %default]', default=6),
    make_option('-y', '--max_ascii_score', action='store', type='string',
                help='SeqPrep option:'+\
                      ' Maximum quality score / ascii code allowed to appear within'+\
                      ' joined pairs output.' 
                      ' [default: %default]', default='J'),
    make_option('-n', '--min_frac_match', action='store', type='float',
                help='SeqPrep option:'+\
                      ' Minimum allowed fraction of matching bases required to join reads'+\
                      ',  [default: %default]', default=0.9),
    make_option('-g', '--max_good_mismatch', action='store', type='float',
                help='SeqPrep option:'+\
                      ' Maximum allowed mis-matched but high quality bases allowed'+\
                      ' to join reads. ' + '[default: %default]', default=0.2),
    make_option('-6', '--phred_64', action='store_true',
                help='SeqPrep option:'+\
                      ' If input reads are in phred+64 format. Output will '\
                      'always be phred+33. [default: %default]',
                        default=False)]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    if not (opts.pe_join_method in join_method_constructors or 
            opts.pe_join_method in join_method_names):
       option_parser.error(\
        'Invalid paired-end join method: %s. \nValid choces are: %s'\
        %opts.pe_join_method,' '.join(join_method_constructors.keys() +
                             join_method_names.keys()))

    forward_reads_fp = opts.forward_reads_fp 
    reverse_reads_fp = opts.reverse_reads_fp
    pe_join_method = opts.pe_join_method
   
    if opts.output_dir: # user specified output directory
        output_dir = abspath(opts.output_dir)
    else: # default output dir to location of infile
        output_dir = dirname(abspath(forward_reads_fp))
    
    create_dir(output_dir, fail_on_exist=False)

    join_func = join_method_names[pe_join_method]

    paths = join_func(forward_reads_fp, reverse_reads_fp, working_dir=output_dir)
   
    if opts.index_reads_fp:
        index_reads = opts.index_reads_fp
        assembly_fp = paths['Assembled']
        
        write_synced_barcodes_fastq(assembly_fp,index_reads)


if __name__ == "__main__":
    main()

