#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight and Micah Hamady"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Micah Hamady", "Greg Caporaso",
               "Kyle Bittinger", "Jesse Stombaugh", "William Walters",
               "Jens Reeder", "Jose Antonio Navas Molina",
               "Jai Ram Rideout"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"


from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from sys import stderr
from qiime.split_libraries import preprocess

options_lookup = get_options_lookup()

# Define the minimum quality score default here so that we can check whether
# the user supplied -s on the command line without also providing -q.
min_qual_score_default = 25

script_info = {}
script_info['brief_description'] = """Split libraries according to barcodes\
 specified in mapping file"""
script_info['script_description'] = """Since newer sequencing technologies\
 provide many reads per run (e.g. the 454 GS FLX Titanium series can produce\
 400-600 million base pairs with 400-500 base pair read lengths) researchers\
 are now finding it useful to combine multiple samples into a single 454 run.\
 This multiplexing is achieved through the application of a\
 pyrosequencing-tailored nucleotide barcode design (described in (Parameswaran\
 et al., 2007)). By assigning individual, unique sample specific barcodes,\
 multiple sequencing runs may be performed in parallel and the resulting reads\
 can later be binned according to sample. The script %prog performs this task,\
 in addition to several quality filtering steps including user defined cut-offs\
 for: sequence lengths; end-trimming; minimum quality score. To summarize, by\
 using the fasta, mapping, and quality files, the program %prog will parse\
 sequences that meet user defined quality thresholds and then rename each read\
 with the appropriate Sample ID, thus formatting the sequence data for\
 downstream analysis. If a combination of different sequencing technologies are\
 used in any particular study, %prog can be used to perform the\
 quality-filtering for each library individually and the output may then be\
 combined.

Sequences from samples that are not found in the mapping file (no corresponding\
 barcode) and sequences without the correct primer sequence will be excluded.\
 Additional scripts can be used to exclude sequences that match a given\
 reference sequence (e.g. the human genome; exclude_seqs_by_blast.py)\
 and/or sequences that are flagged as chimeras (identify_chimeric_seqs.py).
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Standard Example:""",
                                    """Using a single 454 run, which contains a single FASTA, QUAL, and mapping\
 file while using default parameters and outputting the data into the Directory\
 "Split_Library_Output":""",
                                    """%prog -m Mapping_File.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual\
 -o Split_Library_Output/"""))
script_info[
    'script_usage'].append(("""Multiple FASTA and QUAL Files Example:""",
                            """For the case where there are multiple FASTA and QUAL files, the user can run\
 the following comma-separated command as long as there are not duplicate\
 barcodes listed in the mapping file:""",
                            """%prog -m Mapping_File.txt -f 1.TCA.454Reads.fna,2.TCA.454Reads.fna\
 -q 1.TCA.454Reads.qual,2.TCA.454Reads.qual\
 -o Split_Library_Output_comma_separated/"""))
script_info['script_usage'].append(("""Duplicate Barcode Example:""",
                                    """An example of this situation would be a study with 1200 samples. You wish to\
 have 400 samples per run, so you split the analysis into three runs and reuse\
 barcoded primers (you only have 600). After initial analysis you determine a\
 small subset is underrepresented (<500 sequences per samples) and you boost\
 the number of sequences per sample for this subset by running a fourth run.\
 Since the same sample IDs are in more than one run, it is likely that some\
 sequences will be assigned the same unique identifier by %prog when it is run\
 separately on the four different runs, each with their own barcode file. This\
 will cause a problem in file concatenation of the four different runs into a\
 single large file. To avoid this, you can use the '-n' parameter which defines\
 a start index for %prog. From experience, most FLX runs (when combining both\
 files for a single plate) will have 350,000 to 650,000 sequences. Thus, if Run\
 1 for %prog uses '-n 1000000', Run 2 uses '-n 2000000', etc., then you are\
 guaranteed to have unique identifiers after concatenating the results of\
 multiple FLX runs. With newer technologies you will just need to make sure\
 that your start index spacing is greater than the potential number of\
 sequences.

To run %prog, you will need two or more (depending on the number of times the\
 barcodes were reused) separate mapping files (one for each Run, for example\
 one for Run1 and another one for Run2), then you can run %prog using the FASTA\
 and mapping file for Run1 and FASTA and mapping file for Run2. Once you have\
 run split libraries on each file independently, you can concatenate (e.g.\
 using the 'cat' command) the sequence files that were generated by %prog. You\
 can also concatenate the mapping files, since the barcodes are not necessary\
 for downstream analyses, unless the same sample IDs are found in multiple\
 mapping files.

Run %prog on Run 1:""",
                                    """%prog -m Mapping_File.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual\
 -o Split_Library_Run1_Output/ -n 1000000"""))
script_info['script_usage'].append(("""""",
                                    """Run %prog on Run 2. The resulting FASTA files from Run 1 and Run 2 can then\
 be concatenated using the 'cat' command (e.g.\
 cat Split_Library_Run1_Output/seqs.fna Split_Library_Run2_Output/seqs.fna >\
 Combined_seqs.fna) and used in downstream analyses.""",
                                    """%prog -m Mapping_File.txt -f 2.TCA.454Reads.fna -q 2.TCA.454Reads.qual\
 -o Split_Library_Run2_Output/ -n 2000000"""))
script_info['script_usage'].append(("""Barcode Decoding Example:""",
                                    """The standard barcode types supported by %prog are golay (Length: 12 NTs) and\
 hamming (Length: 8 NTs). For situations where the barcodes are of a different\
 length than golay and hamming, the user can define a generic barcode type "-b"\
 as an integer, where the integer is the length of the barcode used in the study.

Note: When analyzing large datasets (>100,000 seqs), users may want to use a\
 generic barcode type, even for length 8 and 12 NTs, since the golay and\
 hamming decoding processes can be computationally intensive, which causes the\
 script to run slow. Barcode correction can be disabled with the -c option if\
 desired.

For the case where the 8 base pair barcodes were used, you can use the\
 following command:""",
                                    """%prog -m Mapping_File_8bp_barcodes.txt -f 1.TCA.454Reads.fna\
  -q 1.TCA.454Reads.qual -o split_Library_output_8bp/ -b 8"""))
script_info['script_usage'].append(("""Linkers and Primers:""",
                                    """The linker and primer sequence (or all the degenerate possibilities) are\
 associated with each barcode from the mapping file. If a barcode cannot be\
 identified, all the possible primers in the mapping file are tested to find a\
 matching sequence. Using truncated forms of the same primer can lead to\
 unexpected results for rare circumstances where the barcode cannot be\
 identified and the sequence following the barcode matches multiple primers.

In many cases, sequence reads are long enough to sequence through the reverse\
 primer and sequencing adapter.  To remove these primers and all following\
 sequences, the -z option can be used.  By default, this option is set to\
 'disable'.  If it is set to 'truncate_only', split_libraries will trim the\
 primer and any sequence following it if the primer is found.  If the\
 'truncate_remove' option is set, %prog will trim the primer if found, and will\
 not write the sequence if the primer is not found. The allowed mismatches for\
 the reverse primer are set with the --reverse_primer_mismatches parameter\
 (default 0).  To use reverse primer removal, one must include a\
 'ReversePrimer' column in the mapping file, with the reverse primer recorded\
 in the 5' to 3' orientation.

Example reverse primer removal, where primers are trimmed if found, and\
 sequence is written unchanged if not found.  Mismatches are increased to 1\
 from the default 0:""",
                                    """%prog -m Mapping_File_reverse_primer.txt -f 1.TCA.454Reads.fna\
 -q 1.TCA.454Reads.qual -o split_libraries_output_revprimer/\
 --reverse_primer_mismatches 1 -z truncate_only"""))

script_info['output_description'] = """Three files are generated by %prog:

1. .fna file (e.g. seqs.fna) - This is a FASTA file containing all sequences\
 which meet the user-defined parameters, where each sequence identifier now\
 contains its corresponding sample id from mapping file.

2. histograms.txt- This contains the counts of sequences with a particular\
 length.

3. split_library_log.txt - This file contains a summary of the %prog analysis.\
 Specifically, this file includes information regarding the number of sequences\
 that pass quality control (number of seqs written) and how these are\
 distributed across the different samples which, through the use of bar-coding\
 technology, would have been pooled into a single 454 run. The number of\
 sequences that pass quality control will depend on length restrictions, number\
 of ambiguous bases, max homopolymer runs, barcode check, etc. All of these\
 parameters are summarized in this file. If raw sequences do not meet the\
 specified quality thresholds they will be omitted from downstream analysis.\
 Since we never see a perfect 454 sequencing run, the number of sequences\
 written should always be less than the number of raw sequences. The number of\
 sequences that are retained for analysis will depend on the quality of the 454\
 run itself in addition to the default data filtering thresholds in the %prog\
 script. The default parameters (minimum quality score = 25, minimum/maximum\
 length = 200/1000, no ambiguous bases allowed, no mismatches allowed in primer\
 sequence) can be adjusted to meet the user's needs.
"""
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
                help='names of qual files, comma-delimited [default: %default]'),
    make_option('-r', '--remove_unassigned', action='store_true',
                help='DEPRECATED: pass --retain_unassigned_reads to keep ' +
                'unassigned reads  [default: %default]'),

    make_option('-l', '--min_seq_length', dest='min_seq_len',
                type='int', default=200,
                help='minimum sequence length, in nucleotides [default: %default]'),

    make_option('-L', '--max_seq_length', dest='max_seq_len',
                type='int', default=1000,
                help='maximum sequence length, in nucleotides [default: %default]'),

    make_option('-t', '--trim_seq_length', dest='trim_seq_len',
                action='store_true',
                help='calculate sequence lengths after trimming primers and barcodes' +
                ' [default: %default]', default=False),

    make_option('-s', '--min_qual_score', type='int', default=None,
                help='min average qual score allowed in read [default: %d]' %
                min_qual_score_default),

    make_option('-k', '--keep_primer', action='store_true',
                help='do not remove primer from sequences', default=False),

    make_option('-B', '--keep_barcode', action='store_true',
                help='do not remove barcode from sequences', default=False),

    make_option('-a', '--max_ambig', type='int', default=6,
                help='maximum number of ambiguous bases [default: %default]'),

    make_option('-H', '--max_homopolymer', type='int', default=6,
                help='maximum length of homopolymer run [default: %default]'),

    make_option('-M', '--max_primer_mismatch', dest='max_primer_mm',
                type='int', default=0,
                help='maximum number of primer mismatches [default: %default]'),

    make_option('-b', '--barcode_type', default='golay_12', type='string',
                help=
                'barcode type, hamming_8, golay_12, variable_length (will ' +
                'disable any barcode correction if variable_length set), or a ' +
                'number representing the length of the barcode, such as -b 4. ' +
                ' [default: %default]'),

    make_option('-o', '--dir_prefix', default='.', type='new_dirpath',
                help='directory prefix for output files [default: %default]'),

    make_option('-e', '--max_barcode_errors', dest='max_bc_errors',
                default=1.5, type='float',
                help='maximum number of errors in barcode [default: %default]'),

    make_option('-n', '--start_numbering_at', dest='start_index',
                default=1, type='int',
                help='seq id to use for the first sequence [default: %default]'),

    make_option("--retain_unassigned_reads", default=False,
                action='store_true', help='retain sequences which are ' +
                'Unassigned in the output sequence file[default: %default]'),

    make_option('-c', '--disable_bc_correction', default=False,
                action='store_true', help='Disable attempts to find nearest ' +
                'corrected barcode.  Can improve performance. [default: %default]'),

    make_option('-w', '--qual_score_window', dest="qual_score_window",
                type='int', default=0,
                action='store', help='Enable sliding window test of quality ' +
                'scores.  If the average score of a continuous set of w nucleotides ' +
                'falls below the threshold (see -s for default), the sequence is ' +
                'discarded. A good value would be 50. 0 (zero) means no filtering. ' +
                'Must pass a .qual file (see -q parameter) if this ' +
                'functionality is enabled.  Default behavior for this function ' +
                'is to truncate the sequence at the beginning of the poor quality ' +
                'window, and test for minimal length (-l parameter) of the resulting ' +
                'sequence. [default: %default]'),

    make_option('-g', '--discard_bad_windows', dest="discard_bad_windows",
                default=False,
                action='store_true', help='If the qual_score_window option (-w) is ' +
                'enabled, this will override the default truncation behavior and ' +
                'discard any sequences where a bad window is found. ' +
                ' [default: %default]'),

    make_option('-p', '--disable_primers', default=False,
                action='store_true', help='Disable primer usage when demultiplexing.' +
                '  Should be enabled for unusual circumstances, such as analyzing ' +
                'Sanger sequence data generated with different primers. ' +
                ' [default: %default]'),

    make_option('-z', '--reverse_primers', default="disable", type='string',
                action='store', help='Enable removal of the reverse primer and ' +
                'any subsequence sequence from the end of each read.  To enable ' +
                'this, there has to be a "ReversePrimer" column in the mapping file. ' +
                "Primers a required to be in IUPAC format and written in the 5' to " +
                " 3' direction.  Valid options are 'disable', 'truncate_only', " +
                "and 'truncate_remove'.  'truncate_only' will remove the primer " +
                "and subsequent sequence data from the output read and will not " +
                "alter output of sequences where the primer cannot be found. " +
                "'truncate_remove' will flag sequences where the primer cannot " +
                "be found to not be written and will record the quantity of such " +
                "failed sequences in the log file. [default: %default]"),

    make_option('--reverse_primer_mismatches', default=0, type='int',
                action='store', help='Set number of allowed mismatches for ' +
                'reverse primers (option -z). [default: %default]'),

    make_option('-d', '--record_qual_scores', default=False,
                action='store_true', help='Enables recording of quality scores for ' +
                'all sequences that are recorded.  If this option is enabled, a file ' +
                'named seqs_filtered.qual will be created in the output directory, ' +
                'and will contain the same sequence IDs in the seqs.fna file and ' +
                'sequence quality scores matching the bases present in the seqs.fna ' +
                'file. [default: %default]'),

    make_option('-i', '--median_length_filtering', default=None, type='string',
                action='store', help='Disables minimum and maximum sequence length ' +
                'filtering, and instead calculates the median sequence length ' +
                'and filters the sequences based upon the number of median absolute ' +
                'deviations specified by this parameter.  Any sequences with lengths ' +
                'outside the number of deviations will be removed. ' +
                '[default: %default]'),

    make_option('-j', '--added_demultiplex_field',
                action='store', default=None, type='string',
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
                ' [default: %default]'),

    make_option('-x', '--truncate_ambi_bases',
                action='store_true', default=False,
                help='Enable to truncate at the first "N" character encountered in ' +
                'the sequences.  This will disable testing for ambiguous bases ' +
                '(-a option) [default: %default]')]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.qual_score_window and not opts.qual_fnames:
        option_parser.error('To enable sliding window quality test (-w), ' +
                            '.qual files must be included.')

    if opts.record_qual_scores and not opts.qual_fnames:
        option_parser.error('To enable recording of truncated quality ' +
                            'scores, one must supply quality score files.')

    min_qual_score = opts.min_qual_score
    if min_qual_score is None:
        min_qual_score = min_qual_score_default
    else:
        if not opts.qual_fnames:
            option_parser.error('To specify a minimum quality score for '
                                'reads, one must supply quality score files.')

    mapping_file = opts.map_fname

    try:
        m = open(mapping_file, "U")
    except IOError:
        raise IOError('Unable to open mapping file %s ' % mapping_file +
                      'Please check filepath and read permissions.')
    fasta_files = set(opts.fasta_fnames)
    if opts.qual_fnames:
        qual_files = set(opts.qual_fnames)
    else:
        qual_files = set()

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

    preprocess(fasta_files, qual_files, mapping_file,
               barcode_type=opts.barcode_type,
               starting_ix=opts.start_index,
               min_seq_len=opts.min_seq_len,
               max_seq_len=opts.max_seq_len,
               min_qual_score=min_qual_score,
               keep_barcode=opts.keep_barcode,
               keep_primer=opts.keep_primer,
               max_ambig=opts.max_ambig,
               max_primer_mm=opts.max_primer_mm,
               trim_seq_len=opts.trim_seq_len,
               dir_prefix=opts.dir_prefix,
               max_bc_errors=opts.max_bc_errors,
               max_homopolymer=opts.max_homopolymer,
               retain_unassigned_reads=opts.retain_unassigned_reads,
               #remove_unassigned = opts.remove_unassigned,
               attempt_bc_correction=not opts.disable_bc_correction,
               qual_score_window=opts.qual_score_window,
               disable_primer_check=opts.disable_primers,
               reverse_primers=opts.reverse_primers,
               reverse_primer_mismatches=opts.reverse_primer_mismatches,
               record_qual_scores=opts.record_qual_scores,
               discard_bad_windows=opts.discard_bad_windows,
               median_length_filtering=opts.median_length_filtering,
               added_demultiplex_field=opts.added_demultiplex_field,
               truncate_ambi_bases=opts.truncate_ambi_bases)

if __name__ == "__main__":
    main()
