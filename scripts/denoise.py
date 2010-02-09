#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

"""Denoising of 454 *.sff.txt files with PyroNoise"""


__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Pre-release"
 
from os import makedirs
from os.path import exists, splitext, split
from optparse import make_option

from cogent.core.alignment import SequenceCollection

from qiime.util import parse_command_line_parameters
from qiime.pyronoise import  pyroNoise_otu_picker

script_description = """Denoise a flowgram file (a .sff.txt, the output of sffinfo)."""

script_usage = """
 Denoise flowgrams in file 454Reads.sff.txt:
   pyronoise.py -i 454Reads.sff.txt 

 Denoise flowgrams in file 454Reads.sff.txt using 2 cores
 on your machine in parallel (requires mpirun):
    pyronoise.py -n 2 -i 454Reads.sff.txt"""

required_options = [\
    make_option('-i','--input_file', action='store',
                type='string', dest='sff_fp',
                help='path to flowgram file (*.sff.txt)')
    ]

optional_options = [\
    make_option('-o','--output_dir', action='store',
                type='string', dest='output_dir',
                help='path to output directory '+
                '[default: %default]',
                default="denoised_seqs/"),
    make_option('-k','--keep_intermediates', action='store_true',
                 dest='keep', default=False,
                 help='Do not delete intermediate files -- '+
                 'useful for debuggingn '+\
                    '[default: %default]'),
    make_option('-c','--cut-off', action='store',
                type='float', dest='cut_off',
                help='cut-off value (passed to pyroNoise) '+\
                    '[default: %default]',
                default = 0.05),
    make_option('-s','--precision', action='store',
                 type='float', dest='precision',
                 help='precision (passed to pyroNoise)'+\
                     '[default: %default]',
                 default=15.0),
    make_option('-n','--num_cpus', action='store',
                type='int', dest='num_cpus',
                help='number of CPUs '+\
                    '[default: %default]',
                default=1)
    ]

def main():
    """run PyroNoise on input flowgrams"""
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    if (not opts.sff_fp or (opts.sff_fp and not exists(opts.sff_fp))):
        option_parser.error(('Flowgram file path does not exist:\n %s \n'+\
                                 'Pass a valid one via -i.')% opts.sff_fp)

    input_seqs_dir, input_seqs_filename = split(opts.sff_fp)
    #split off .txt
    input_seqs_basename, ext = splitext(input_seqs_filename)
    #split off .sff
    input_seqs_basename, ext = splitext(input_seqs_basename)

    outdir = opts.output_dir

    if (not exists(outdir)):
        try:
            makedirs(outdir)
        except OSError:
            #re-raise error, but slightly more informative 
            raise OSError,"Could not create output directory "+outdir
    log_fh=None
    if (opts.verbose):
        try:
            log_fh = open(outdir+"/pyronoise.log", "w")
        except IOError:
            raise IOError,"Could not open log file: %s" % (outdir+"pyronoise.log")

        #write params to log file
        # should have a general framework for this in util...
        log_fh.write("Input file: %s\n"% opts.sff_fp)
        log_fh.write("output path: %s\n"% outdir)

    centroids, cluster_mapping = pyroNoise_otu_picker(open(opts.sff_fp, "U"),
                                                      outdir, opts.num_cpus, log_fh, opts.keep,
                                                      opts.precision, opts.cut_off)

    # store mapping file and centroids
    result_otu_path = '%s/%s_otus.txt' % (outdir, input_seqs_basename)
    of = open(result_otu_path,'w')
    for i,cluster in cluster_mapping.iteritems():
        of.write('%s\t%s\n' % (i,'\t'.join(cluster)))
    of.close()
    
    result_fasta_path = '%s/%s.fasta' % (outdir, input_seqs_basename)
    of = open(result_fasta_path,'w')
    of.write(SequenceCollection(centroids).toFasta()+"\n")

if __name__ == "__main__":
    main()
