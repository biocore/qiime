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
    "Emily TerAvest"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"

from string import upper
from itertools import izip
from os.path import join
from os import rename
from collections import defaultdict
from operator import itemgetter
from gzip import GzipFile

from numpy import array
from skbio.parse.sequences import parse_fasta

from qiime.parse import QiimeParseError, MinimalQualParser
from qiime.hamming import decode_barcode_8
from qiime.golay import decode as decode_golay_12
from qiime.check_id_map import process_id_map
from qiime.barcode import correct_barcode

""" This library contains the code for demultiplexing 454 data.  Apart from
    barcode correction/mismatch counts, this code does not quality
    filter sequences.

    Output includes a seqs.fna file in the specified output directory with
    the sequences demultiplexed according to barcode or other demultiplexing
    fields.  If quality score file(s) are supplied, these are also written
    with the same labels, and barcode sequence removed.  A log with details
    about sequences assigned per sample and unassigned counts will be created
    as well.  If write_unassigned_reads is enabled, all sequences that can
    not be demultiplexed will be written to seqs_not_assigned.fna (and .qual
    if qual score file(s) are supplied)."""


def process_files_and_demultiplex_sequences(mapping_file,
                                            fasta_files,
                                            qual_files,
                                            output_dir="./",
                                            keep_barcode=False,
                                            barcode_type='golay_12',
                                            max_bc_errors=0.5,
                                            start_index=1,
                                            write_unassigned_reads=False,
                                            disable_bc_correction=False,
                                            added_demultiplex_field=None,
                                            save_barcode_frequencies=False):
    """ Handles file IO, calls main demultiplexing function

    mapping_file:  filepath to metadata mapping file.
    fasta_files:  list of fasta filepaths
    qual_files: list of qual filepaths
    output_dir:  output directory to write all output files.
    keep_barcode:  If True, will not remove barcode from output files.
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    max_bc_errors:  Number of changes allowed for error correcting barcodes,
     for generic barcodes, specifies the number of mismatches allowed.
    start_index:  Specifies the first number used to enumerate output sequences.
    write_unassigned_reads:  If True, will write sequences that could not be
     demultiplexed into a separate output file.
    disable_bc_correction:  Only tests for exact matches to barcodes.
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.
    save_barcode_frequencies:  Saves the frequencies of barcode sequences in
     a separate output file.
    """

    file_data = {}

    fasta_files = [get_infile(fasta_f) for fasta_f in fasta_files]
    qual_files = [get_infile(qual_f) for qual_f in qual_files]

    file_data['fasta_files'] = fasta_files
    file_data['qual_files'] = qual_files
    file_data['mapping_file'] = open(mapping_file, "U")

    file_data['demultiplexed_seqs_f'] = open(join(output_dir,
                                                  "demultiplexed_seqs.fna.incomplete"), "w")
    if qual_files:
        file_data['demultiplexed_qual_f'] = open(join(output_dir,
                                                      "demultiplexed_seqs.qual.incomplete"), "w")
    if write_unassigned_reads:
        file_data['unassigned_seqs_f'] = open(join(output_dir,
                                                   "unassigned_seqs.fna.incomplete"), "w")
        if qual_files:
            file_data['unassigned_qual_f'] =\
                open(join(output_dir, "unassigned_seqs.qual.incomplete"), "w")

    log_data, bc_freqs, seq_counts, corrected_bc_count =\
        demultiplex_sequences(file_data, keep_barcode, barcode_type,
                              max_bc_errors, start_index, write_unassigned_reads,
                              disable_bc_correction, added_demultiplex_field)

    final_log_data = process_log_data(log_data, seq_counts, mapping_file,
                                      fasta_files, qual_files, corrected_bc_count, keep_barcode, barcode_type,
                                      max_bc_errors, start_index, write_unassigned_reads, disable_bc_correction,
                                      added_demultiplex_field)

    log_file = open(join(output_dir, "demultiplex_fasta.log"), "w")
    log_file.write("\n".join(final_log_data))

    if save_barcode_frequencies:
        bcs_sorted_list = process_bc_freqs(bc_freqs)
        bc_freqs_f = open(join(output_dir, "barcode_freqs.txt"), "w")
        bc_freqs_f.write("Barcode frequencies\n")
        bc_freqs_f.write("\n".join(bcs_sorted_list))

    # Rename .incomplete files to .fna/.qual files

    rename(file_data['demultiplexed_seqs_f'].name, join(output_dir,
                                                        "demultiplexed_seqs.fna"))
    if qual_files:
        rename(file_data['demultiplexed_qual_f'].name, join(output_dir,
                                                            "demultiplexed_seqs.qual"))
    if write_unassigned_reads:
        rename(file_data['unassigned_seqs_f'].name, join(output_dir,
                                                         "unassigned_seqs.fna"))
        if qual_files:
            rename(file_data['unassigned_qual_f'].name, join(output_dir,
                                                             "unassigned_seqs.qual"))


def demultiplex_sequences(file_data,
                          keep_barcode=False,
                          barcode_type='golay_12',
                          max_bc_errors=1.5,
                          start_index=1,
                          write_unassigned_reads=False,
                          disable_bc_correction=False,
                          added_demultiplex_field=None):
    """ Main program function for demultiplexing fasta sequence data

    file_data:  dict of open file objects, contains input fasta, qual, and
     mapping files, and output filepaths for partially demultiplexed fasta
     and qual files, and unassigned sequence output file.
    keep_barcode:  If True, will not remove barcode from output files.
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    max_bc_errors:  Number of changes allowed for error correcting barcodes,
     for generic barcodes, specifies the number of mismatches allowed.
    start_index:  Specifies the first number used to enumerate output sequences.
    write_unassigned_reads:  If True, will write sequences that could not be
     demultiplexed into a separate output file.
    disable_bc_correction:  Only tests for exact matches to barcodes.
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.
    """

    header, mapping_data = check_map(file_data['mapping_file'], barcode_type,
                                     added_demultiplex_field)

    ids_bcs_added_field = get_ids_bcs_added_field(header, mapping_data,
                                                  barcode_type, added_demultiplex_field)

    bc_lens = get_bc_lens(ids_bcs_added_field)

    # Get all bcs in a list here to save computation later
    all_bcs = [curr_bc[0] for curr_bc in ids_bcs_added_field.keys()]

    log_data, bc_freqs, seq_counts, corrected_bc_count =\
        assign_seqs(
            file_data, ids_bcs_added_field, bc_lens, all_bcs, keep_barcode,
            barcode_type, max_bc_errors, start_index, write_unassigned_reads,
            disable_bc_correction, added_demultiplex_field)

    return log_data, bc_freqs, seq_counts, corrected_bc_count


def assign_seqs(file_data,
                ids_bcs_added_field,
                bc_lens,
                all_bcs,
                keep_barcode=False,
                barcode_type="golay_12",
                max_bc_errors=1.5,
                start_index=1,
                write_unassigned_reads=False,
                disable_bc_correction=False,
                added_demultiplex_field=None):
    """ Demultiplexes, writes seqs/qual files, returns log data

    file_data:  dict of open file objects, contains input fasta, qual, and
     mapping files, and output filepaths for partially demultiplexed fasta
     and qual files, and unassigned sequence output file.
    ids_bcs_added_field: dict of (barcode,added_demultiplex): SampleID
    bc_lens:  Lengths of all barcodes from largest to smallest.
    all_bcs:  List of all barcode sequences.
    keep_barcode:  If True, will not remove barcode from output files.
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    max_bc_errors:  Number of changes allowed for error correcting barcodes,
     for generic barcodes, specifies the number of mismatches allowed.
    start_index:  Specifies the first number used to enumerate output sequences.
    write_unassigned_reads:  If True, will write sequences that could not be
     demultiplexed into a separate output file.
    disable_bc_correction:  Only tests for exact matches to barcodes.
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.
    save_barcode_frequencies:  Saves the frequencies of barcode sequences in
     a separate output file.
    """

    log_data = initialize_log_data(ids_bcs_added_field)
    bc_freqs = defaultdict(int)

    seq_counts = 0
    enum_val = start_index
    corrected_bc_count = [0, 0]

    if file_data['qual_files']:
        for curr_fasta, curr_qual in zip(file_data['fasta_files'],
                                         file_data['qual_files']):
            for fasta_data, qual_data in izip(parse_fasta(curr_fasta),
                                              MinimalQualParser(curr_qual, full_header=True)):

                seq_counts += 1
                fasta_label, fasta_seq = fasta_data
                qual_label, qual_seq = qual_data

                bc, corrected_bc, num_errors, added_field =\
                    get_demultiplex_data(ids_bcs_added_field,
                                         fasta_label, fasta_seq, bc_lens, all_bcs, barcode_type,
                                         max_bc_errors, disable_bc_correction, added_demultiplex_field)

                bc_freqs[bc] += 1

                sample_id, log_id, bc_corrected_result =\
                    get_output_ids(ids_bcs_added_field,
                                   corrected_bc, num_errors, added_field, max_bc_errors,
                                   enum_val)
                if bc_corrected_result == 'corrected':
                    corrected_bc_count[0] += 1
                if bc_corrected_result == 'not_corrected':
                    corrected_bc_count[1] += 1

                label_line = get_label_line(sample_id, fasta_label, bc,
                                            corrected_bc, num_errors)

                if sample_id.startswith("Unassigned") and\
                        write_unassigned_reads:
                    write_fasta_line(file_data['unassigned_seqs_f'],
                                     fasta_seq, label_line, True, len(bc))
                    write_qual_line(file_data['unassigned_qual_f'],
                                    list(qual_seq), label_line, True, len(bc))
                elif not sample_id.startswith("Unassigned"):
                    write_fasta_line(file_data['demultiplexed_seqs_f'],
                                     fasta_seq, label_line, keep_barcode, len(bc))
                    write_qual_line(file_data['demultiplexed_qual_f'],
                                    list(qual_seq), label_line, keep_barcode, len(bc))

                if log_id:
                    log_data[log_id] += 1

                enum_val += 1

    else:
        for curr_fasta in file_data['fasta_files']:
            for fasta_label, fasta_seq in parse_fasta(curr_fasta):
                seq_counts += 1
                bc, corrected_bc, num_errors, added_field =\
                    get_demultiplex_data(ids_bcs_added_field,
                                         fasta_label, fasta_seq, bc_lens, all_bcs, barcode_type,
                                         max_bc_errors, disable_bc_correction, added_demultiplex_field)

                bc_freqs[bc] += 1

                sample_id, log_id, bc_corrected_result =\
                    get_output_ids(ids_bcs_added_field,
                                   corrected_bc, num_errors, added_field, max_bc_errors,
                                   enum_val)

                if bc_corrected_result == 'corrected':
                    corrected_bc_count[0] += 1
                if bc_corrected_result == 'not_corrected':
                    corrected_bc_count[1] += 1

                label_line = get_label_line(sample_id, fasta_label, bc,
                                            corrected_bc, num_errors)

                if sample_id.startswith("Unassigned") and\
                        write_unassigned_reads:
                    write_fasta_line(file_data['unassigned_seqs_f'],
                                     fasta_seq, label_line, True, len(bc))
                elif not sample_id.startswith("Unassigned"):
                    write_fasta_line(file_data['demultiplexed_seqs_f'],
                                     fasta_seq, label_line, keep_barcode, len(bc))

                if log_id:
                    log_data[log_id] += 1

                enum_val += 1

    return log_data, bc_freqs, seq_counts, corrected_bc_count


def get_output_ids(ids_bcs_added_field,
                   corrected_bc,
                   num_errors,
                   added_field,
                   max_bc_errors=1.5,
                   enum_val=1):
    """ Returns SampleID to write to output fasta/qual files

    ids_bcs_added_field: dict of (barcode,added_demultiplex): SampleID
    corrected_bc: Corrected barcode sequence, None if correction not possible
    num_errors: Barcode errors/mismatches for current barcode
    added_field:  Current added demultiplex field, None if not found or if
     no added demultiplex fields used.
    max_bc_errors: Max number of barcode errors/mismatches allowed
    enum_val: Value to follow SampleID and underscore
    """

    bc_corrected_flag = None

    if added_field is None:
        curr_added_field = ''
    else:
        curr_added_field = added_field
    if corrected_bc is None:
        curr_bc = ''
    else:
        curr_bc = corrected_bc
    log_id = ""
    if num_errors > max_bc_errors:
        sample_id = "Unassigned_%d" % enum_val
        bc_corrected_flag = 'not_corrected'
    else:
        try:
            base_sample_id =\
                ids_bcs_added_field[(curr_bc, curr_added_field)]
            sample_id = "%s_%d" % (base_sample_id, enum_val)
            if corrected_bc:
                log_id += "%s" % corrected_bc
            if corrected_bc and added_field:
                log_id += ","
            if added_field:
                log_id += "%s" % added_field
            # Correction for case of only having a single sample with no
            # demultiplexing data
            if log_id:
                log_id += ","
            log_id += "%s" % base_sample_id

            if num_errors > 0:
                bc_corrected_flag = 'corrected'
        except KeyError:
            sample_id = "Unassigned_%d" % enum_val

    # Returns sample_id for writing fasta label, log_id for output counts,
    # and bc_corrected_flag to count number of barcodes that are corrected
    # bc_corrected_flag is None for no corrections, 'not_corrected' if errors
    # exceed max allowed, and 'corrected' if errors are within allowed limit.
    # None can either mean that no barcode correction was necessary due to
    # an exact match or no correction was possible at all.
    return sample_id, log_id, bc_corrected_flag


def initialize_log_data(ids_bcs_added_field):
    """ Initializes log data, so that zero count Samples can be recorded

    ids_bcs_added_field: dict of (barcode,added_demultiplex): SampleID
    """

    log_data = {}

    for curr_key in ids_bcs_added_field.keys():
        base_key = ""
        if curr_key[0]:
            base_key += curr_key[0] + ","
        if curr_key[1]:
            base_key += curr_key[1] + ","
        base_key += ids_bcs_added_field[curr_key]
        log_data[base_key] = 0

    return log_data


def get_label_line(sample_id,
                   fasta_label,
                   bc,
                   corrected_bc,
                   num_errors):
    """ Returns line to use for fasta/qual output label

    sample_id:  enumerated SampleID
    fasta_label:  Original fasta label
    bc: Original barcode sequence
    corrected_barcode:  Corrected barcode sequence
    num_errors:  Number of errors/mismatches in barcode sequence
    """

    orig_label = fasta_label.split()[0]

    final_label = "%s %s orig_bc=%s new_bc=%s bc_diffs=%d" %\
        (sample_id, orig_label, bc, corrected_bc, num_errors)

    return final_label


def write_fasta_line(demultiplexed_seqs_f,
                     fasta_seq,
                     label_line,
                     keep_barcode,
                     bc_len):
    """ Writes fasta label, seq to demultiplexed_seqs_f file

    demultiplexed_seqs_f:  open file object to write label/seq to.
    fasta_seq:  current fasta sequence
    label_line:  fasta label to write
    keep_barcode:  If True, do not slice out barcode from sequence
    bc_len:  Length of barcode, used to slice out from sequence
    """

    if keep_barcode:
        final_seq = fasta_seq
    else:
        final_seq = fasta_seq[bc_len:]

    demultiplexed_seqs_f.write(">%s\n" % label_line)
    demultiplexed_seqs_f.write("%s\n" % final_seq)


def write_qual_line(demultiplexed_qual_f,
                    qual_seq,
                    label_line,
                    keep_barcode,
                    bc_len):
    """ Writes quality score sequence out in proper format

    demultiplexed_qual_f:  open file object to write label/qual scores to.
    qual_seq:  current qual sequence, list of scores.
    label_line:  qual label to write
    keep_barcode:  If True, do not slice out barcode from sequence
    bc_len:  Length of barcode, used to slice out from sequence
    """

    if keep_barcode:
        final_seq = qual_seq
    else:
        final_seq = qual_seq[bc_len:]

    qual_line_size = 60

    current_qual_scores_lines = []
    # Quality score format is a string of 60 base calls, followed by a
    # newline, until the last N bases are written
    for slice in range(0, len(final_seq), qual_line_size):

        current_segment = final_seq[slice:slice + qual_line_size]
        current_qual_scores_lines.append(" ".join(map(str, current_segment)))

    demultiplexed_qual_f.write(">%s\n" % label_line)
    demultiplexed_qual_f.write('\n'.join(current_qual_scores_lines))
    demultiplexed_qual_f.write('\n')


def get_demultiplex_data(ids_bcs_added_field,
                         fasta_label,
                         fasta_seq,
                         bc_lens,
                         all_bcs,
                         barcode_type="golay_12",
                         max_bc_errors=1.5,
                         disable_bc_correction=False,
                         added_demultiplex_field=None):
    """ Attempts to find bc in a given sequence and added demultiplex field

    ids_bcs_added_field:  dict of (barcode,added_demultiplex): SampleID
    fasta_label:  full fasta label, needed for added demultiplex field
    fasta_seq:  current fasta sequence
    bc_lens: size of barcodes, from largest to smallest
    all_bcs: List of all barcode sequences.
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    max_bc_errors:  Number of changes allowed for error correcting barcodes,
     for generic barcodes, specifies the number of mismatches allowed.
    disable_bc_correction:  Only tests for exact matches to barcodes.
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.
    """

    # To allow for variable length barcodes, need to step down from largest
    # possible barcode size to smallest, break on an exact match, otherwise
    # return None, 0, added label indices if matched.
    for bc_len in bc_lens:
        curr_bc = fasta_seq[0:bc_len]
        corrected_bc, num_errors, added_field = get_curr_bc_added_field(
            curr_bc,
            ids_bcs_added_field, fasta_label, all_bcs, barcode_type,
            disable_bc_correction, added_demultiplex_field)

        # Escape if exact hit found for barcode, only matters for variable
        # length barcodes.  Need special case for variable length barcodes
        # that have overlapping sequences and added demultiplex field.
        if added_field:
            if (corrected_bc, added_field) in ids_bcs_added_field.keys():
                break
        elif corrected_bc is not None:
            break

    return curr_bc, corrected_bc, num_errors, added_field


def get_curr_bc_added_field(curr_bc,
                            ids_bcs_added_field,
                            fasta_label,
                            all_bcs,
                            barcode_type="golay_12",
                            disable_bc_correction=False,
                            added_demultiplex_field=None):
    """ Attempts to correct barcode, get added demultiplex data

    curr_bc: current barcode sequence to attempt correction with
    ids_bcs_added_field: dict of (barcode,added_demultiplex): SampleID
    fasta_label:  Current full fasta label
    all_bcs:  List of all barcode sequences.
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    disable_bc_correction:  Only tests for exact matches to barcodes.
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.
    """

    if added_demultiplex_field:
        added_field = get_added_demultiplex_field(ids_bcs_added_field,
                                                  fasta_label, added_demultiplex_field)
    else:
        added_field = None

    if disable_bc_correction:
        num_errors = 0
        corrected_bc = get_exact_bc_matches(curr_bc, all_bcs)
    else:
        corrected_bc, num_errors = attempt_bc_correction(curr_bc,
                                                         all_bcs, barcode_type)

    return corrected_bc, num_errors, added_field


def attempt_bc_correction(curr_bc,
                          all_bcs,
                          barcode_type="golay_12"):
    """ Gets corrected barcode and number of errors

    curr_bc: current barcode sequence to attempt correction with
    all_bcs: List of all barcode sequences.
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    """

    # First check for exact matches
    corrected_bc = get_exact_bc_matches(curr_bc, all_bcs)
    if corrected_bc:
        return corrected_bc, 0

    if barcode_type == "golay_12":
        corrected_bc, num_errors = decode_golay_12(curr_bc)
    elif barcode_type == "hamming_8":
        corrected_bc, num_errors = decode_barcode_8(curr_bc)
    elif barcode_type == 0:
        corrected_bc, num_errors = ('', 0)
    else:
        corrected_bc, num_errors = correct_barcode(curr_bc, all_bcs)

    return corrected_bc, num_errors


def get_exact_bc_matches(curr_bc,
                         all_bcs):
    """ Checks existing barcodes for an exact match, returns None if not found

    curr_bc: current barcode sequence to attempt correction with
    all_bcs: List of all barcode sequences
    """

    if curr_bc in all_bcs:
        return curr_bc
    else:
        return None


def get_added_demultiplex_field(ids_bcs_added_field,
                                fasta_label,
                                added_demultiplex_field):
    """ Returns matches to added demultiplex field from fasta label

    ids_bcs_added_field: dict of (barcode,added_demultiplex): SampleID
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.
    """

    for curr_bc_added_field in ids_bcs_added_field.keys():
        curr_added_field = curr_bc_added_field[1]
        if added_demultiplex_field == 'run_prefix':
            curr_label_slice = fasta_label[0:len(curr_added_field)]
            if curr_label_slice == curr_added_field:
                return curr_label_slice
        else:
            # Skip if current added field not in fasta label
            if curr_added_field not in fasta_label:
                continue
            curr_label_slice =\
                fasta_label.split(
                    added_demultiplex_field + "=")[1][0:len(curr_added_field)]
            if curr_label_slice == curr_added_field:
                return curr_label_slice

    # If nothing found, return None
    return None


def check_map(mapping_file,
              barcode_type="golay_12",
              added_demultiplex_field=None):
    """ Gets header, mapping data, halts execution if there are errors

    mapping_file:  list of lines of metadata mapping file
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.
    """

    if barcode_type == 0:
        has_barcodes = False
        var_len_barcodes = False
    elif barcode_type == 'variable_length':
        has_barcodes = True
        var_len_barcodes = True
    else:
        has_barcodes = True
        var_len_barcodes = False

    header, mapping_data, run_description, errors, warnings = \
        process_id_map(mapping_file, has_barcodes=has_barcodes,
                       disable_primer_check=True,
                       added_demultiplex_field=added_demultiplex_field,
                       variable_len_barcodes=var_len_barcodes)

    # Need to specifically detect varied length barcodes, otherwise won't know
    # how much of sequence to slice off for barcode reads
    for warning in warnings:
        if "differs than length" in warning:
            raise ValueError("Detected variable length barcodes, if these " +
                             "are being used, use -b variable_length")
    # Halt on errors, as these are serious problems with mapping file.
    # These include non-DNA characters in the barcodes, duplicate
    # barcodes or duplicate barcodes/added demultiplex fields, duplicate
    # SampleIDs, or header problems.
    if errors:
        raise ValueError("Errors found in mapping file, please check " +
                         "mapping file with validate_mapping_file.py")

    return header, mapping_data


def get_ids_bcs_added_field(header,
                            mapping_data,
                            barcode_type="golay_12",
                            added_demultiplex_field=None):
    """ Creates dict of data for SampleID, barcodes, and added field.

    header:  list of strings of header data.
    mapping_data:  list of lists of mapping data
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.

    Dict is of format key = (barcode,added_demultiplex): SampleID
    """
    sample_id_ix = header.index("SampleID")
    bc_ix = header.index("BarcodeSequence")
    if added_demultiplex_field:
        added_demultiplex_ix = header.index(added_demultiplex_field)

    ids_bcs_added_field = {}

    for line in mapping_data:

        if barcode_type == 0:
            curr_bc = ''
        else:
            curr_bc = line[bc_ix]
        if added_demultiplex_field:
            curr_added_field = line[added_demultiplex_ix]
        else:
            curr_added_field = ''

        ids_bcs_added_field[(upper(curr_bc), curr_added_field)] =\
            line[sample_id_ix]

    return ids_bcs_added_field


def get_bc_lens(ids_bcs_added_field):
    """ Get list of barcode lens, sorted from largest to smallest

    ids_bcs_added_field: dict of (barcode,added_demultiplex): SampleID
    """

    bc_index = 0
    bc_lens = [len(curr_bc[bc_index])
               for curr_bc in ids_bcs_added_field.keys()]
    bc_lens = list(set(bc_lens))
    bc_lens.sort(reverse=True)
    return bc_lens


def get_infile(filename):
    """Returns filehandle, allowing gzip input."""
    if filename.endswith(".gz"):
        fin = GzipFile(filename, "rb")
    else:
        fin = open(filename, "U")
    return fin


def process_log_data(log_data,
                     seq_counts,
                     mapping_file,
                     fasta_files,
                     qual_files,
                     corrected_bc_count,
                     keep_barcode=False,
                     barcode_type="golay_12",
                     max_bc_errors=1.5,
                     start_index=1,
                     write_unassigned_reads=False,
                     disable_bc_correction=False,
                     added_demultiplex_field=None,
                     save_barcode_frequencies=False):
    """ Formats log data in list of lines for writing

    log_data: dict of barcode,added_demultiplex: sequence counts
    seq_counts:  Total counts of sequences demultiplexed
    mapping_file:  filepath to metadata mapping file.
    fasta_files:  list of open fasta file objects
    qual_files: list of open qual file objects
    corrected_bc_count: two element list, first element is counts of corrected
     barcodes, the second is counts of barcodes with too many errors to correct
    keep_barcode:  If True, will not remove barcode from output files.
    barcode_type:  Specified barcode, can be golay_12, hamming_8,
     variable_length, or an integer specifying length.
    max_bc_errors:  Number of changes allowed for error correcting barcodes,
     for generic barcodes, specifies the number of mismatches allowed.
    start_index:  Specifies the first number used to enumerate output sequences.
    write_unassigned_reads:  If True, will write sequences that could not be
     demultiplexed into a separate output file.
    disable_bc_correction:  Only tests for exact matches to barcodes.
    added_demultiplex_field:  Uses data supplied in metadata mapping field
     and demultiplexes according to data in fasta labels.
    save_barcode_frequencies:  Saves the frequencies of barcode sequences in
     a separate output file.
    """

    final_log_data = ["demultiplex_fasta.py log data\n"]
    final_log_data.append("Metadata mapping file:\t%s" % mapping_file)
    final_log_data.append("Input FASTA file(s):\t%s" %
                          ",".join([curr_file.name for curr_file in fasta_files]))
    if qual_files:
        final_log_data.append("Input QUAL file(s):\t%s" %
                              ",".join([curr_file.name for curr_file in qual_files]))
    final_log_data.append("Total sequences in input files:\t%d" % seq_counts)
    final_log_data.append("Retain barcode:\t%s" % keep_barcode)
    final_log_data.append("Barcode type:\t%s" % barcode_type)
    final_log_data.append("Max barcode error/mismatches allowed:\t%s" %
                          max_bc_errors)
    final_log_data.append("Starting sequence identifier:\t%d" % start_index)
    final_log_data.append(
        "Write unassigned reads:\t%s" %
        write_unassigned_reads)
    final_log_data.append("Disable barcode correction:\t%s" %
                          disable_bc_correction)
    final_log_data.append("Added demultiplex field:\t%s" %
                          added_demultiplex_field)
    final_log_data.append("Save barcode frequencies:\t%s\n" %
                          save_barcode_frequencies)

    final_log_data.append("Barcodes corrected/not corrected:\t%d/%d" %
                          (corrected_bc_count[0], corrected_bc_count[1]))

    # Make a list of SampleID/Counts/Barcodes/Added demultiplex for sorting
    id_counts_bcs = []
    counts_col = []
    for curr_key in log_data.keys():
        id_counts_bcs.append((curr_key.split(',')[-1], log_data[curr_key],
                              ','.join(curr_key.split(',')[0:-1])))
        counts_col.append(int(log_data[curr_key]))

    counts_col = array(counts_col)

    final_log_data.append("Number of samples in mapping file:\t%d" %
                          len(counts_col))
    final_log_data.append("Sample count min/max/mean:\t%d / %d / %3.2f" %
                          (counts_col.min(), counts_col.max(), counts_col.mean()))

    id_counts_bcs = sorted(id_counts_bcs, key=itemgetter(1), reverse=True)

    final_log_data.append("Sample\tSequence Count\tBarcode/Added Demultiplex")
    for curr_id in id_counts_bcs:
        final_log_data.append("%s\t%s\t%s" %
                              (curr_id[0], curr_id[1], curr_id[2]))

    final_log_data.append("Seqs written\t%d" % counts_col.sum())
    final_log_data.append("Percent of input seqs written\t%3.2f" %
                          (counts_col.sum() / seq_counts))

    return final_log_data


def process_bc_freqs(bc_freqs):
    """ Sorts barcode frequencies according to counts, returns list

    bc_freqs:  dictionary of barcodes:counts
    """

    bcs_list = []
    for curr_key in bc_freqs.keys():
        bcs_list.append((curr_key, int(bc_freqs[curr_key])))

    bcs_list = sorted(bcs_list, key=itemgetter(1), reverse=True)

    sorted_bcs = []
    for curr_bc in bcs_list:
        sorted_bcs.append("%s\t%d" % (curr_bc[0], curr_bc[1]))

    return sorted_bcs
