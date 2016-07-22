#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = [
    "William Walters",
    "Rob Knight",
    "Micah Hamady",
    "Greg Caporaso",
    "Kyle Bittinger",
    "Jesse Stombaugh",
    "Jens Reeder",
    "Jose Antonio Navas Molina"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option, create_dir

from qiime.demultiplex_fasta import process_files_and_demultiplex_sequences

options_lookup = get_options_lookup()
script_info = {}
script_info[
    'brief_description'] = """Demultiplex fasta data according to barcode sequences or data supplied in fasta labels."""
script_info[
    'script_description'] = """Using barcodes and/or data from fasta labels provided in a mapping file, will demultiplex sequences from an input fasta file.  Barcodes will be removed from the sequences in the output fasta file by default.  If a quality scores file is supplied, the quality score file will be truncated to match the output fasta file.  The default barcode type are 12 base pair Golay codes.  Alternative barcodes allowed are 8 base pair Hamming codes, variable_length, or generic barcodes of a specified length.  Generic barcodes utilize mismatch counts for correction.  One can also use an added demultiplex field (-j option) to specify data in the fasta labels that can be used alone or in conjunction with barcode sequences for demultiplexing.  All barcode correction is disabled when variable length barcodes are used."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Standard Example:""",
     """Using a single 454 run, which contains a single FASTA, QUAL, and mapping file while using default parameters and outputting the data into the Directory "demultiplexed_output":""",
     """%prog -m Mapping_File_golay.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual -o demultiplexed_output/"""))
script_info['script_usage'].append(
    ("""""",
     """For the case where there are multiple FASTA and QUAL files, the user can run the following command as long as there are not duplicate barcodes listed in the mapping file:""",
     """%prog -m Mapping_File_golay.txt -f 1.TCA.454Reads.fna,2.TCA.454Reads.fna -q 1.TCA.454Reads.qual,2.TCA.454Reads.qual -o demultiplexed_output_comma_separated/"""))
script_info['script_usage'].append(("""Duplicate Barcode Example:""", """An example of this situation would be a study with 1200 samples. You wish to have 400 samples per run, so you split the analysis into three runs with and reuse barcodes (you only have 600). After initial analysis you determine a small subset is underrepresented (<500 sequences per samples) and you boost the number of sequences per sample for this subset by running a fourth run. Since the same sample IDs are in more than one run, it is likely that some sequences will be assigned the same unique identifier by %prog when it is run separately on the four different runs, each with their own barcode file. This will cause a problem in file concatenation of the four different runs into a single large file. To avoid this, you can use the '-n' parameter which defines a start index for %prog fasta label enumeration. From experience, most 454 runs (when combining both files for a single plate) will have 350,000 to 650,000 sequences. Thus, if Run 1 for %prog uses '-n 1000000', Run 2 uses '-n 2000000', etc., then you are guaranteed to have unique identifiers after concatenating the results of multiple 454 runs. With newer technologies you will just need to make sure that your start index spacing is greater than the potential number of sequences.

To run %prog, you will need two or more (depending on the number of times the barcodes were reused) separate mapping files (one for each Run, for example one Run1 and another one for Run2), then you can run %prog using the FASTA and mapping file for Run1 and FASTA and mapping file for Run2. Once you have independently run demultiplex_fasta on each file, followed by quality filtering, you can concatenate (cat) the sequence files generated. You can also concatenate the mapping files, since the barcodes are not necessary for downstream analyses, unless the same sample ids are found in multiple mapping files.

Run %prog on Run 1:""", """%prog -m Mapping_File1.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual -o demultiplexed_output_Run1/ -n 1000000"""))
script_info['script_usage'].append(
    ("""""",
     """Run demultiplex_fasta on Run 2:""",
     """%prog -m Mapping_File2.txt -f 2.TCA.454Reads.fna -q 2.TCA.454Reads.qual -o demultiplexed_output_Run2/ -n 2000000"""))
script_info['script_usage'].append(("""Barcode Decoding Example:""", """The standard barcode types supported by %prog are golay (Length: 12 NTs) and hamming (Length: 8 NTs). For situations where the barcodes are of a different length than golay and hamming, the user can define a generic barcode type "-b" as an integer, where the integer is the length of the barcode used in the study.

For the case where the generic 8 base pair barcodes were used, you can use the following command:""", """%prog -m Mapping_File_8bp_barcodes.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual -o demultiplexed_output_8bp_barcodes/ -b 8"""))
script_info['script_usage'].append(
    ("""""",
     """To use the run prefix at the beginning of the fasta label for demultiplexing, there has to be a field in the mapping file labeled "run_prefix", and can be used by the following command:""",
     """%prog -m Mapping_File_run_prefix.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual -o demultiplexed_output_run_prefix/ -j run_prefix"""))

script_info['output_description'] = """Four files can be generated by %prog

1. seqs.fna - This contains the fasta sequences, demultiplexed according to barcodes and/or added demultiplexed field.

2. demultiplexed_sequences.log - Contains details about demultiplexing stats

3. seqs.qual - If quality score file(s) are supplied, these will be truncated to match the seqs.fna file after barcode removal if such is enabled.

4. seqs_not_assigned.fna - If --retain_unassigned_reads is enabled, will write all sequences that can not be demultiplexed to this file.  Also will create a seqs_not_assigned.qual file if quality file supplied."""

script_info['required_options'] = [
    make_option('-m', '--map', dest='map_fname', type='existing_filepath',
                help='name of mapping file. NOTE: Must contain a header' +
                ' line indicating SampleID in the first column and' +
                ' BarcodeSequence in the second,' +
                ' LinkerPrimerSequence in the third.'),
    make_option(
        '-f', '--fasta', dest='fasta_fnames', type='existing_filepaths',
        help='names of fasta files, comma-delimited')]

script_info['optional_options'] = [
    make_option('-q', '--qual', dest='qual_fnames', type='existing_filepaths',
                help='file paths of qual files, comma-delimited [default: %default]'),

    make_option('-B', '--keep_barcode', action='store_true',
                help='do not remove barcode from sequences', default=False),

    make_option('-b', '--barcode_type', default='golay_12',
                help=
                'barcode type, hamming_8, golay_12, variable_length (will ' +
                'disable any barcode correction if variable_length set), or a ' +
                'number representing the length of the barcode, such as -b 4. ' +
                'The max barcode errors (-e) should be lowered for short barcodes. ' +
                ' [default: %default]'),

    make_option('-o', '--dir_prefix', default='.', type='new_dirpath',
                help='directory prefix for output files [default: %default]'),

    make_option('-e', '--max_barcode_errors', dest='max_bc_errors',
                default=1.5, type='float',
                help='maximum number of errors in barcode.  If using generic barcodes' +
                ' every 0.5 specified counts as a primer mismatch.' +
                ' [default: %default]'),

    make_option('-n', '--start_numbering_at', dest='start_index',
                default=1, type='int',
                help='seq id to use for the first sequence [default: %default]'),

    make_option("--retain_unassigned_reads", default=False,
                action='store_true', help='retain sequences which can not be ' +
                'demultiplexed in a seperate output sequence file [default: %default]'),

    make_option('-c', '--disable_bc_correction', default=False,
                action='store_true', help='Disable attempts to find nearest ' +
                'corrected barcode.  Can improve performance. [default: %default]'),

    make_option('-F', '--save_barcode_frequencies', default=False,
                action='store_true', help='Save frequences of barcodes as they ' +
                'appear in the given sequences.  Sorts in order of largest to ' +
                'smallest.  Will do nothing if barcode type is 0 or variable_length. ' +
                ' [default: %default]'),

    make_option('-j', '--added_demultiplex_field',
                action='store', default=None,
                help='Use -j to add a field to use in the mapping file as an ' +
                'additional demultiplexing option to the barcode.  All combinations ' +
                'of barcodes and the values in these fields must be unique. The ' +
                'fields must contain values that can be parsed from the fasta labels ' +
                'such as "plate=R_2008_12_09".  In this case, "plate" would be the ' +
                'column header and "R_2008_12_09" would be the field data (minus ' +
                'quotes) in the mapping file.  To use the run prefix from the fasta ' +
                'label, such as ">FLP3FBN01ELBSX", where "FLP3FBN01" is generated ' +
                'from the run ID, use "-j run_prefix" and set the run prefix to ' +
                'be used as the data under the column headerr "run_prefix". ' +
                ' [default: %default]')]


script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    mapping_file = opts.map_fname
    fasta_files = opts.fasta_fnames
    output_dir = opts.dir_prefix
    keep_barcode = opts.keep_barcode
    barcode_type = opts.barcode_type
    max_bc_errors = opts.max_bc_errors
    start_index = opts.start_index
    write_unassigned_reads = opts.retain_unassigned_reads
    disable_bc_correction = opts.disable_bc_correction
    added_demultiplex_field = opts.added_demultiplex_field
    save_barcode_frequencies = opts.save_barcode_frequencies

    # Test filepaths
    try:
        m = open(mapping_file, "U")
    except IOError:
        raise IOError('Unable to open mapping file %s ' % mapping_file +
                      'Please check filepath and read permissions.')

    if opts.qual_fnames:
        qual_files = opts.qual_fnames
    else:
        qual_files = []

    for q in qual_files:
        try:
            test_qual_file = open(q, "U")
            test_qual_file.close()
        except IOError:
            raise IOError('Unable to open file %s ' % q + ' please check ' +
                          'filepath and read permissions.')

    for f in fasta_files:
        try:
            test_fasta_file = open(f, "U")
            test_fasta_file.close()
        except IOError:
            raise IOError('Unable to open file %s ' % f + ' please check ' +
                          'filepath and read permissions.')

    # Test for duplicates in fasta and qual files
    if len(set(fasta_files)) != len(fasta_files):
        raise ValueError("Duplicate fasta files found, please check file " +
                         "file paths specified")
    if len(set(qual_files)) != len(qual_files):
        raise ValueError("Duplicate qual files found, please check file " +
                         "file paths specified")
    # Should have same number of fasta and qual files
    if qual_files:
        if len(fasta_files) != len(qual_files):
            raise ValueError("Number of fasta and qual files differ, please " +
                             "check files specified with -f and -q.")

    # Check barcode type for valid selection
    valid_bc_types = ['golay_12', 'hamming_8', 'variable_length']
    if barcode_type not in valid_bc_types:
        try:
            barcode_type = int(barcode_type)
            if barcode_type < 0:
                raise ValueError
        except ValueError:
            raise ValueError('Barcode types must be golay_12, hamming_8, ' +
                             'variable_length, or a positive integer value representing the ' +
                             'length of the barcode.')
    # Attempts at barcode correction disabled for variable length barcodes
    if barcode_type == 'variable_length':
        disable_bc_correction = True
    if barcode_type == 'variable_length' or barcode_type == 0:
        save_barcode_frequencies = False

    create_dir(output_dir)

    process_files_and_demultiplex_sequences(mapping_file, fasta_files,
                                            qual_files, output_dir, keep_barcode, barcode_type, max_bc_errors,
                                            start_index, write_unassigned_reads, disable_bc_correction,
                                            added_demultiplex_field, save_barcode_frequencies)

if __name__ == "__main__":
    main()
