#!/usr/bin/env python
# File created Feb 1 2012

__author__ = "William Anton Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Anton Walters"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Anton Walters"
__email__ = "william.a.walters@gmail.com"

from qiime.util import parse_command_line_parameters, get_options_lookup,\
    make_option, create_dir
from qiime.validate_demultiplexed_fasta import validate_fasta

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Checks a fasta file to verify if it has  been properly demultiplexed, i.e., it is in QIIME compatible format."""
script_info[
    'script_description'] = """Checks file is a valid fasta file, does not contain gaps ('.' or '-' characters), contains only valid nucleotide characters, no fasta label is duplicated, SampleIDs match those in a provided mapping file, fasta labels are formatted to have SampleID_X as normally generated by QIIME demultiplexing, and the BarcodeSequence/LinkerPrimerSequences are not found in the fasta sequences.  Optionally this script can also verify that the SampleIDs in the fasta sequences are also present in the tip IDs of a provided newick tree file, can test for equal sequence lengths across all sequences, and can test that all SampleIDs in the mapping file are represented in the fasta file labels."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example:""",
     """ """,
     """ validate_demultiplexed_fasta.py -f seqs.fasta -m Mapping_File.txt"""))
script_info['output_description'] = """"""
script_info['required_options'] = [
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='name of mapping file. NOTE: Must contain a header' +
                ' line indicating SampleID in the first column and' +
                ' BarcodeSequence in the second,' +
                ' LinkerPrimerSequence in the third.  If no barcode or ' +
                ' linkerprimer sequence is present, leave data fields empty.'),

    options_lookup['fasta_as_primary_input']
]
script_info['optional_options'] = [
    make_option('-o', '--output_dir', default='.', type='new_dirpath',
                help='directory prefix for output files [default: %default]'),
    make_option('-t', '--tree_fp', default=None, type='existing_filepath',
                help='path to the tree file; ' +
                'Needed to test if sequence IDs are a subset or exact match to the ' +
                'tree tips, options -s and -e  [default: %default]'),
    make_option('-s', '--tree_subset', default=False, action='store_true',
                help='Determine if sequence IDs are a subset of the tree tips, ' +
                'newick tree must be passed with the -t option. [default: %default]'),
    make_option('-e', '--tree_exact_match', default=False, action='store_true',
                help='Determine if sequence IDs are an exact match to tree tips, ' +
                'newick tree must be passed with the -t option. [default: %default]'),
    make_option('-l', '--same_seq_lens', default=False, action='store_true',
                help='Determine if sequences are all the same length. ' +
                '[default: %default]'),
    make_option('-a', '--all_ids_found', default=False, action='store_true',
                help='Determine if all SampleIDs provided in the mapping file ' +
                'are represented in the fasta file labels. [default: %default]'),
    make_option('-b', '--suppress_barcode_checks', default=False,
                action='store_true', help='Suppress barcode checks ' +
                '[default: %default]'),
    make_option('-p', '--suppress_primer_checks', default=False,
                action='store_true', help='Suppress primer checks ' +
                '[default: %default]')

]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose

    input_fasta_fp = opts.input_fasta_fp
    mapping_fp = opts.mapping_fp
    output_dir = opts.output_dir
    tree_fp = opts.tree_fp
    tree_subset = opts.tree_subset
    tree_exact_match = opts.tree_exact_match
    same_seq_lens = opts.same_seq_lens
    all_ids_found = opts.all_ids_found

    create_dir(output_dir)

    # Test optional filepaths and requirements
    try:
        test_mapping_fp = open(mapping_fp, "U")
        test_mapping_fp.close()
    except IOError:
        raise IOError("Unable to open mapping file, please check "
                      "filepath and read permissions.")

    if tree_fp:
        try:
            test_tree_fp = open(tree_fp, "U")
            test_tree_fp.close()
        except IOError:
            raise IOError("Unable to open provided tree filepath, please " +
                          "filepath and permissions.")

    if tree_subset or tree_exact_match:
        if not tree_fp:
            raise ValueError('Must provide tree filepath if -s or -e options ' +
                             'are enabled.')

    validate_fasta(
        input_fasta_fp, mapping_fp, output_dir, tree_fp, tree_subset,
        tree_exact_match, same_seq_lens, all_ids_found,
        opts.suppress_barcode_checks, opts.suppress_primer_checks)


if __name__ == "__main__":
    main()
