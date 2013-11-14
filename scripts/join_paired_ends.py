#!/usr/bin/env python
# file: join_paired_ends.py
# Using make_phylogeny as a guide.

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
                                    join_method_constructors)
from qiime.util import (parse_command_line_parameters, get_options_lookup, 
                        make_option, load_qiime_config, create_dir)
import gzip

options_lookup = get_options_lookup()
qiime_config = load_qiime_config()

script_info={}
script_info['brief_description']= """Joins paired-end Illumina reads."""
script_info['script_description'] = """This script take forward and reverse Illumina reads and joins them based on the method chosen. Currently, there are four methods that can be selected by the user:

1. fastq-join - Erik Aronesty, 2011. ea-utils : "Command-line tools for processing biological sequencing data" (http://code.google.com/p/ea-utils)

2. SeqPrep - (https://github.com/jstjohn/SeqPrep)

3. FLASh - Magoc & Salzberg (2011) Bioinformatics (http://ccb.jhu.edu/software/FLASH/)

4. PandaSeq - Masella et al. (2012) PANDAseq: paired-end assembler for illumina sequences. BMC Bioinformatics 13:31. (https://github.com/neufeld/pandaseq)

"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Join paired-ends with \'fastq-join\':""","""This is the default method to join paired-end Illumina data:""","""%prog -f $PWD/forward_reads.fastq -r reverse_reads.fastq"""))
script_info['script_usage'].append(("""Join paired-ends with \'SeqPrep\':""","""Produces similar output to the \'fastq-join\' but returns gzipped data by default (this option can be changed).""","""%prog -m SeqPrep -f $PWD/forward_reads.fastq -r reverse_reads.fastq"""))
script_info['script_usage'].append(("""Join paired-ends with \'FLASh\':""","""This method is multi-thread capable but should only be used if the paired-ends have already been quality trimmed, or are not highly overlapping.""","""%prog -m flash -f $PWD/forward_reads.fastq -r reverse_reads.fastq"""))
script_info['script_usage'].append(("""Join paired-ends with \'PandaSeq\':""","""This method is also multi-thread capable. Unlike the previous methods, the output from \'PandaSeq\' can not be used for downstream quality filtering. \'PandaSeq\' changes the meaning of the quality scores in regions of overlap. This may cause unexcpected behaviour with downstream quality filters.""","""%prog -m pandaseq -f $PWD/forward_reads.fastq -r reverse_reads.fastq"""))
script_info['output_description'] = """All paired-end joining software will return a joined / merged / assembled paired-end fastq file. Depending on the method chosen, additional files may be written to the user-specified output directory. 

The following will always be returned:
\"..._joined.fastq\" - This is a FASTQ file containing the joined paired-end sequences.

Additional output given the method chosen via this script:
1. fastq-join, SeqPrep, & FLASh will output: 
   \"..._un1.fastq\" - unassembled / unjoined reads1 output
   \"..._un2.fastq\" - unassembled / unjoined reads2 output

"""
script_info['required_options'] = [\
    make_option('-f','--forward_reads_fp',type="existing_filepath",
                 dest='forward_reads_fp',
                 help='Path to read input forward reads in FASTQ format.'),
    make_option('-r','--reverse_reads_fp',type="existing_filepath",
                 dest='reverse_reads_fp',
                 help='Path to read input reverse reads in FASTQ format.')]
    #make_option('-b','--index_reads_fp',type="existing_filepath",
    #            help='Path to read the barcode / index reads in FASTQ format.'+\
    #                ' Will be filtered based on surviving joined pairs.')
script_info['optional_options'] = [\
    make_option('-m', '--pe_join_method', action='store', type='choice',
                choices=list(join_method_names.keys()),
                help='Method to use for joining paired-ends. Valid choices'+\
                      ' are: ' + ', '.join(join_method_names.keys())+\
                      ' [default: %default]', default='fastq-join'),
    make_option('-o', '--output_dir', action='store', type='new_dirpath',\
                help='Path to store '+\
                      'result file [default: <JOINED_METHOD>_joined]'),
     make_option('-t', '--threads', action='store', type='new_dirpath',\
                help='Number of cpus to use. '+\
                      'Only applicable when the method used is'+\
                      'FLASh\' or \'PandaSeq\' ')]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    if not (opts.join_method in join_method_constructors or 
            opts.join_method in join_method_names):
       option_parser.error(\
        'Invalid paired-end join method: %s. \nValid choces are: %s'\
        %opts.pe_join_method,' '.join(join_method_constructors.keys() +
                             join_method_names.keys()))

    forward_reads_fp = opts.forward_reads_fp 
    reverse_reads_fp = opts.reverse_reads_fp
    pe_join_method = opts.pe_join_method
    
    if opt.output_dir:
        output_dir = opt.output_dir
    else:
        fpath,ext = splitext(forward_reads_fp)[0]
        output_dir = fpath + pe_join_method + '_joined'
    
    create_dir(output_dir, fail_on_exists=False)

    join_func = join_method_names(pe_join_method)

    paths = join_func(forward_reads_fp, reverse_reads_fp, working_dir=output_dir)

if __name__ == "__main__":
    main()

