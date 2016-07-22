#!/usr/bin/env python
# File created February 29, 2012
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from os.path import join, basename

from skbio.parse.sequences import parse_fasta
from skbio.sequence import DNA

from qiime.split_libraries import local_align_primer_seq
from qiime.check_id_map import process_id_map


def get_rev_primer_seqs(mapping_fp):
    """ Parses mapping file to get dictionary of SampleID:Rev primer
    mapping_fp:  mapping filepath
    """
    hds, mapping_data, run_description, errors, warnings = \
        process_id_map(mapping_fp, has_barcodes=False,
                       disable_primer_check=True)

    if errors:
        for curr_err in errors:
            if curr_err.startswith("Duplicate SampleID"):
                raise ValueError('Errors were found with mapping file, ' +
                                 'please run validate_mapping_file.py to ' +
                                 'identify problems.')

    # create dict of dicts with SampleID:{each header:mapping data}

    id_map = {}

    for curr_data in mapping_data:
        id_map[curr_data[0]] = {}

    for header in range(len(hds)):
        for curr_data in mapping_data:
            id_map[curr_data[0]][hds[header]] = curr_data[header]

    reverse_primers = {}

    for curr_id in id_map.keys():
        try:
            reverse_primers[curr_id] =\
                [str(DNA(curr_rev_primer).rc()) for curr_rev_primer in
                 id_map[curr_id]['ReversePrimer'].split(',')]
        except KeyError:
            raise KeyError("Reverse primer not found in mapping file, " +
                           "please include a 'ReversePrimer' column.")

    # Check for valid reverse primers
    # Will have been detected as warnings from mapping file
    for curr_err in errors:
        if curr_err.startswith("Invalid DNA sequence detected"):
            raise ValueError("Problems found with reverse primers, please " +
                             "check mapping file with validate_mapping_file.py")

    return reverse_primers


def get_output_filepaths(output_dir,
                         fasta_fp):
    """ Returns output fasta filepath and log filepath

    fasta_fp: fasta filepath
    output_dir: output directory
    """
    fasta_extensions = ['.fa', '.fasta', '.fna']

    curr_fasta_out = basename(fasta_fp)
    for fasta_extension in fasta_extensions:
        curr_fasta_out = curr_fasta_out.replace(fasta_extension, '')

    curr_fasta_out += "_rev_primer_truncated.fna"
    output_fp = join(output_dir, curr_fasta_out)
    log_fp = join(output_dir, "rev_primer_truncation.log")

    return output_fp, log_fp


def truncate_rev_primers(fasta_f,
                         output_fp,
                         reverse_primers,
                         truncate_option='truncate_only',
                         primer_mismatches=2):
    """ Locally aligns reverse primers, trucates or removes seqs

    fasta_f:  open file of fasta file
    output_fp: open filepath to write truncated fasta to
    reverse_primers: dictionary of SampleID:reverse primer sequence
    truncate_option: either truncate_only, truncate_remove
    primer_mismatches: number of allowed primer mismatches
    """

    log_data = {
        'sample_id_not_found': 0,
        'reverse_primer_not_found': 0,
        'total_seqs': 0,
        'seqs_written': 0
    }

    for label, seq in parse_fasta(fasta_f):
        curr_label = label.split('_')[0]

        log_data['total_seqs'] += 1

        # Check fasta label for valid SampleID, if not found, just write seq
        try:
            curr_rev_primer = reverse_primers[curr_label]
        except KeyError:
            log_data['sample_id_not_found'] += 1
            output_fp.write('>%s\n%s\n' % (label, seq))
            log_data['seqs_written'] += 1
            continue

        mm_tests = {}
        for rev_primer in curr_rev_primer:

            rev_primer_mm, rev_primer_index =\
                local_align_primer_seq(rev_primer, seq)

            mm_tests[rev_primer_mm] = rev_primer_index

        rev_primer_mm = min(mm_tests.keys())
        rev_primer_index = mm_tests[rev_primer_mm]

        if rev_primer_mm > primer_mismatches:
            if truncate_option == "truncate_remove":
                log_data['reverse_primer_not_found'] += 1
            else:
                log_data['reverse_primer_not_found'] += 1
                log_data['seqs_written'] += 1
                output_fp.write('>%s\n%s\n' % (label, seq))
        else:
            # Check for zero seq length after truncation, will not write seq
            if rev_primer_index > 0:
                log_data['seqs_written'] += 1
                output_fp.write('>%s\n%s\n' % (label, seq[0:rev_primer_index]))

    return log_data


def write_log_file(log_data,
                   log_f):
    """ Writes log file

    log_data: dictionary of details about reverse primer removal
    log_f: open filepath to write log details
    """

    log_f.write("Details for removal of reverse primers\n")
    log_f.write("Original fasta filepath: %s\n" % log_data['fasta_fp'])
    log_f.write("Total seqs in fasta: %d\n" % log_data['total_seqs'])
    log_f.write("Mapping filepath: %s\n" % log_data['mapping_fp'])
    log_f.write("Truncation option: %s\n" % log_data['truncate_option'])
    log_f.write("Mismatches allowed: %d\n" % log_data['primer_mismatches'])
    log_f.write("Total seqs written: %d\n" % log_data['seqs_written'])
    log_f.write("SampleIDs not found: %d\n" % log_data['sample_id_not_found'])
    log_f.write("Reverse primers not found: %d\n" %
                log_data['reverse_primer_not_found'])


def truncate_reverse_primer(fasta_fp,
                            mapping_fp,
                            output_dir=".",
                            truncate_option='truncate_only',
                            primer_mismatches=2):
    """ Main program function for finding, removing reverse primer seqs

    fasta_fp: fasta filepath
    mapping_fp: mapping filepath
    output_dir: output directory
    truncate_option: truncation option, either truncate_only, truncate_remove
    primer_mismatches: Number is mismatches allowed in reverse primer"""

    reverse_primers = get_rev_primer_seqs(open(mapping_fp, "U"))

    output_fp, log_fp = get_output_filepaths(output_dir, fasta_fp)

    log_data = truncate_rev_primers(open(fasta_fp, "U"),
                                    open(
                                        output_fp, "w"), reverse_primers, truncate_option,
                                    primer_mismatches)

    log_data['fasta_fp'] = fasta_fp
    log_data['mapping_fp'] = mapping_fp
    log_data['truncate_option'] = truncate_option
    log_data['primer_mismatches'] = primer_mismatches

    write_log_file(log_data, open(log_fp, "w"))
