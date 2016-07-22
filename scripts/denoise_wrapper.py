#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

"""Denoising of 454 *.sff.txt files"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder", "Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os.path import exists, splitext, split
from qiime.util import make_option
from numpy import array

from burrito.util import ApplicationError

from qiime.util import parse_command_line_parameters, create_dir
from qiime.denoise_wrapper import fast_denoiser
from qiime.parse import parse_mapping_file
from qiime.format import write_Fasta_from_name_seq_pairs

script_info = {}
script_info['brief_description'] = """Denoise a flowgram file"""
script_info[
    'script_description'] = """This script will denoise a flowgram file in .sff.txt format, which is the output of sffinfo."""

script_info['script_usage'] = [
    ("""Example:""",
     """Denoise flowgrams in file 454Reads.sff.txt, discard flowgrams not in seqs.fna, and extract primer from map.txt:""",
     """%prog -i 454Reads.sff.txt -f seqs.fna -m map.txt"""),

    ("""Multi-core Example:""",
     """Denoise flowgrams in file 454Reads.sff.txt using 2 cores on your machine in parallel:""",
     """%prog -n 2 -i 454Reads.sff.txt -f seqs.fna -m map.txt""")
]

script_info[
    'output_description'] = """This script results in a OTU like mapping file along with a sequence file of denoised (FASTA-format). Note that the sequences coming from denoising are no real OTUs, and have to be sent to pick_otus.py if the users wishes to have a defined similarity threshold."""

script_info['required_options'] = [
    make_option('-i', '--input_file', action='store',
                type='existing_filepaths', dest='sff_fps',
                help='path to flowgram files (.sff.txt), ' +
                'comma separated'),

    make_option('-f', '--fasta_file', action='store',
                type='existing_filepath', dest='fasta_fp',
                help='path to fasta file from split_libraries.py')
]

script_info['optional_options'] = [
    make_option('-o', '--output_dir', action='store',
                type='new_dirpath', dest='output_dir',
                help='path to output directory ' +
                '[default: %default]',
                default="denoised_seqs/"),

    make_option('-n', '--num_cpus', action='store',
                type='int', dest='num_cpus',
                help='number of CPUs ' +
                '[default: %default]',
                default=1),

    make_option('--force_overwrite', action='store_true',
                dest='force', default=False,
                help='Overwrite files in output directory ' +
                '[default: %default]'),

    make_option('-m', '--map_fname', action='store',
                type='existing_filepath', dest='map_fname',
                help='name of mapping file, Has to contain ' +
                'field LinkerPrimerSequence. ' +
                '[REQUIRED unless --primer specified]'),

    make_option('-p', '--primer', action='store',
                type='string', dest='primer',
                help='primer sequence ' +
                '[REQUIRED unless --map_fname specified]',
                default=None),

    make_option('--titanium', action='store_true',
                dest='titanium', default=False,
                help='Select Titanium defaults for denoiser, '
                + 'otherwise use FLX defaults ' +
                '[default: %default]')
]

script_info['version'] = __version__


def main():
    """run denoiser on input flowgrams"""
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    sff_files = opts.sff_fps

    for f in sff_files:
        if (not exists(f)):
            option_parser.error(('Flowgram file path does not exist:\n %s \n' +
                                 'Pass a valid one via -i.') % f)
    outdir = opts.output_dir

    create_dir(outdir, fail_on_exist=not opts.force)

    log_fh = None

    if (not (opts.primer or opts.map_fname)):
        raise ApplicationError("Either mapping file or primer required")
    # Read primer from Meta data file if not set on command line
    if not opts.primer:
        mapping_data, header, comments = \
            parse_mapping_file(open(opts.map_fname, "U"))

        index = header.index("LinkerPrimerSequence")
        all_primers = set(array(mapping_data)[:, index])

        if len(all_primers) != 1:
            raise ValueError("Currently only data sets with one primer are allowed.\n" +
                             "Make separate mapping files with only one primer, re-run split_libraries and\n"
                             + "denoise with each split_library output separately.")
        primer = list(all_primers)[0]
        last_char = primer[-1]
        if(last_char not in "ACGT"):
            raise ValueError("We currently do not support primer with " +
                             "degenerate bases at it's 3' end.")

    else:
        primer = opts.primer

    centroids, cluster_mapping = fast_denoiser(opts.sff_fps, opts.fasta_fp,
                                               outdir, opts.num_cpus, primer,
                                               titanium=opts.titanium)

    # store mapping file and centroids
    result_otu_path = '%s/denoised_clusters.txt' % outdir
    of = open(result_otu_path, 'w')
    for i, cluster in cluster_mapping.iteritems():
        of.write('%s\t%s\n' % (str(i), '\t'.join(cluster)))
    of.close()

    result_fasta_path = '%s/denoised_seqs.fasta' % outdir
    oh = open(result_fasta_path, 'w')
    write_Fasta_from_name_seq_pairs(centroids, oh)

if __name__ == "__main__":
    main()
