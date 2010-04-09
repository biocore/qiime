#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

"""Denoising of 454 *.sff.txt files"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Pre-release"
 
from os.path import exists, splitext, split
from optparse import make_option
from numpy import array

from cogent.core.alignment import SequenceCollection

from qiime.util import parse_command_line_parameters, create_dir,\
    handle_error_codes
from qiime.pyronoise import  pyroNoise_otu_picker, fast_denoiser
from qiime.parse import parse_mapping_file

script_info={}
script_info['brief_description']="""Denoise a flowgram file"""
script_info['script_description']="""This script will denoise a flowgram file in  .sff.txt format, which is the output of sffinfo."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Denoise flowgrams in file 454Reads.sff.txt:""","""denoise.py -i 454Reads.sff.txt"""))
script_info['script_usage'].append(("""Multi-core Example:""","""Denoise flowgrams in file 454Reads.sff.txt using 2 cores on your machine in parallel (requires mpirun):""","""denoise.py -n 2 -i 454Reads.sff.txt"""))
script_info['output_description']="""This script results in a OTU mapping file along with a sequence file of denoised (FASTA-format). Note that the sequences coming from denoising are no real OTUs, and have to be sent to pick_otus.py if the users wishes to have a defined similarity threshold. """

script_info['required_options'] = [\
    make_option('-i','--input_file', action='store',
                type='string', dest='sff_fp',
                help='path to flowgram file (.sff.txt)'),

    make_option('-f','--fasta_file', action='store',
                type='string', dest='fasta_fp',
                help='path to fasta file from split_libraries.py')
    ]

script_info['optional_options'] = [\
    make_option('-o','--output_dir', action='store',
                type='string', dest='output_dir',
                help='path to output directory '+
                '[default: %default]',
                default="denoised_seqs/"),

    make_option('--method', action='store',
                type='string', dest='method',
                help='Method to use for denoising. '+
                'Choice of pyronoise or fast'+
                ' [default: %default]',
                default="fast"),

    make_option('-k','--keep_intermediates', action='store_true',
                 dest='keep', default=False,
                 help='Do not delete intermediate files -- '+
                 'useful for debugging '+\
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
                default=1),
    make_option('--force_overwrite', action='store_true',
                 dest='force', default=False,
                 help='Overwrite files in output directory '+\
                    '[default: %default]'),

    make_option('-m','--map_fname', action='store',
                type='string', dest='map_fname',
                help='name of mapping file, Has to contain '+\
                    'field LinkerPrimerSequence.'+\
                    ' [REQUIRED] when method is fast'),
    
    make_option('-p', '--primer',action='store',\
                    type='string',dest='primer',\
                    help='primer sequence '+\
                    '[default: %default]',
                default=None)
    ]

script_info['version'] = __version__

def main():
    """run PyroNoise on input flowgrams"""
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if (not opts.sff_fp or (opts.sff_fp and not exists(opts.sff_fp))):
        option_parser.error(('Flowgram file path does not exist:\n %s \n'+\
                                 'Pass a valid one via -i.')% opts.sff_fp)

    input_seqs_dir, input_seqs_filename = split(opts.sff_fp)
    #split off .txt
    input_seqs_basename, ext = splitext(input_seqs_filename)
    #split off .sff
    input_seqs_basename, ext = splitext(input_seqs_basename)
    outdir = opts.output_dir

    ret_val = create_dir(outdir, handle_errors_externally=True)  
    if ret_val==1:  #dir exists
        if opts.force:
            #do nothing, just overwrite content
            pass
        else:
            print "Directory exists. Use --force to overwrite."
            exit()
    else:
        handle_error_codes(outdir, error_code=ret_val)

    log_fh=None
    if (opts.verbose):
        try:
            log_fh = open(outdir+"/pyronoise.log", "w")
        except IOError:
            raise IOError,"Could not open log file: %s" % (outdir+"/pyronoise.log")

        #write params to log file
        # should have a general framework for this in util...
        log_fh.write("Input file: %s\n"% opts.sff_fp)
        log_fh.write("output path: %s\n"% outdir)

    if opts.method=="pyronoise":
        centroids, cluster_mapping = pyroNoise_otu_picker(open(opts.sff_fp, "U"),
                                                          outdir, opts.num_cpus, log_fh, opts.keep,
                                                          opts.precision, opts.cut_off)
    else:
        #Read primer from Meta data file if not set on command line
        if not opts.primer:
          mapping_data, header, comments = \
              parse_mapping_file(open(opts.map_fname,"U"))
            
          index = header.index("LinkerPrimerSequence")
          all_primers = set(array(mapping_data)[:,index])
          
          if len(all_primers)!= 1:
                raise ValueError,"Currently only data sets with one primer are allowed.\n"+\
                    "Make separate mapping files with only one primer, re-run split_libraries and\n"\
                    +"denoise with each split_library output separately."
          primer = list(all_primers)[0]
          last_char = primer[-1]
          if(last_char not in "ACGT"):
              raise ValueError,"We currently do not support primer with "+\
                  "degenerate bases at it's 3' end."

        else:
            primer=opts.primer

        centroids, cluster_mapping = fast_denoiser(opts.sff_fp,opts.fasta_fp,
                                                   outdir, opts.num_cpus, primer)

    # store mapping file and centroids
    result_otu_path = '%s/%s_denoised_clusters.txt' % (outdir, input_seqs_basename)
    of = open(result_otu_path,'w')
    for i,cluster in cluster_mapping.iteritems():
        of.write('%s\t%s\n' % (str(i),'\t'.join(cluster)))
    of.close()
    
    result_fasta_path = '%s/denoised_seqs.fasta' % outdir
    of = open(result_fasta_path,'w')
    of.write(SequenceCollection(dict(centroids)).toFasta()+"\n")

if __name__ == "__main__":
    main()
