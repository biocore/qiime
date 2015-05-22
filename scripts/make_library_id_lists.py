#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Doug Wendel"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Doug Wendel"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"


from qiime.make_library_id_lists import get_first_id, get_ids
from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from os.path import exists, join
from os import makedirs

script_info = {}
script_info['brief_description'] = """Make library id lists"""
script_info[
    'script_description'] = """Makes a list of the ids corresponding to each library represented in the input fasta file. Assumes that the libraries are the output of split_libraries.py and that they contain the 454 read id for each sequence as is standard in the split_libraries.py output. Produces a separate file for each library."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example:""",
     """Create a list containing library ids for a fasta file (seqs.fna):""",
     """make_library_id_lists.py -i seqs.fna -o results/"""))
script_info[
    'output_description'] = """This script produces a separate file for each library."""
script_info['required_options'] = [
    make_option(
        "-i", "--input_fasta", dest='in_fasta', default=None, type='existing_filepath',
        help="The path to a FASTA file containing input sequences")
]

script_info['optional_options'] = [
    make_option("-s", "--screened_rep_seqs", dest="screened_rep_seqs",
                default=None, type='existing_filepath',
                help="The path to a FASTA file containing screened representative seqs" +
                "[DEFAULT: %default]"),
    make_option("-u", "--otus", dest="otus",
                default=None, type='existing_filepath',
                help="The path to an OTU file mapping OTUs onto rep seqs" +
                "[DEFAULT: %default]"),
    make_option("-o", "--outdir", dest='outdir',
                default='.', type='new_dirpath',
                help=""" The base directory to save results (one file per library)."""),
    make_option("-f", "--field", dest="field", type='int',
                default=1,
                help="Index of space-delimited field to read id from [DEFAULT: %default]"),
    make_option("--debug", dest="debug", action="store_true",
                default=False, help="Show debug output.")
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    options, args = option_parser.parse_args()
    if options.debug:
        print "PRODUCING DEBUG OUTPUT"

    bad_seq_ids = set()
    bad_otu_ids = None

    # if we got a file to screen against, find the relevant ids and delete them
    if options.screened_rep_seqs:
        bad_otu_ids = get_first_id(open(options.screened_rep_seqs, 'U'))
        if not options.otus:
            raise RuntimeError(
                "Must specify an OTU file if performing a screen.")
        for line in open(options.otus, 'U'):
            fields = line.split()
            if fields[0] in bad_otu_ids:
                bad_seq_ids.update(fields[1:])

    if options.debug:
        if bad_otu_ids is not None:
            print "Found %s bad otu ids: %s" % (len(bad_otu_ids), bad_otu_ids)
        print "Found %s bad seq ids: %s" % (len(bad_seq_ids), bad_seq_ids)

    ids = get_ids(open(options.in_fasta, 'U'), options.field, bad_seq_ids,
                  options.debug)

    # add empty unassigned ids for file creation
    if 'Unassigned' not in ids:
        ids['Unassigned'] = []

    if not exists(options.outdir):
        makedirs(options.outdir)
    for k, idlist in ids.items():
        outfile = open(join(options.outdir, k + '.txt'), 'w')
        outfile.write('\n'.join(sorted(idlist)))
        outfile.close()

if __name__ == "__main__":
    main()
