#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@gmail.com"

import numpy as np

from string import upper
from itertools import izip, cycle
from os.path import join
from os import rename
from re import compile

from skbio.parse.sequences import parse_fastq
from skbio.sequence import DNA
from skbio.format.sequences import format_fastq_record

from qiime.check_id_map import process_id_map
from qiime.split_libraries_fastq import (check_header_match_pre180,
                                         check_header_match_180_or_later)
from qiime.parse import is_casava_v180_or_later
from qiime.pycogent_backports.fastq import FastqParseError


def extract_barcodes(fastq1,
                     fastq2=None,
                     output_dir=".",
                     input_type="barcode_single_end",
                     bc1_len=6,
                     bc2_len=6,
                     rev_comp_bc1=False,
                     rev_comp_bc2=False,
                     char_delineator=":",
                     switch_bc_order=False,
                     map_fp=None,
                     attempt_read_orientation=False,
                     disable_header_match=False):
    """ Main program function for extracting barcodes from reads

    fastq1: Open fastq file 1.
    fastq2: None or open fastq file 2.
    output_dir: Directory to write output parses sequences to.
    input_type: Specifies the type of parsing to be done.
    bc1_len: Length of barcode 1 to be parsed from fastq1
    bc2_len: Length of barcode 2 to be parsed from fastq2, or from end of a
     stitched read.
    rev_comp_bc1: If True, reverse complement bc1 before writing.
    rev_comp_bc2: If True, reverse complement bc2 before writing.
    char_delineator: Specify character that immediately precedes the barcode
        for input_type of barcode_in_label.
    switch_bc_order: Normally, barcode 1 will be written first, followed by
        barcode 2 in a combined output fastq file. If True, the order will be
        reversed. Only applies to stitched reads processing, as other barcode
        orders are dictated by the the parameter chosen for the fastq files.
    map_fp: open file object of mapping file, requires a LinkerPrimerSequence
        and ReversePrimer field to be present. Used for orienting reads.
    attempt_read_orientation: If True, will attempt to orient the reads
        according to the forward primers in the mapping file. If primer is
        detected in current orientation, leave the read as is, but if reverse
        complement is detected (or ReversePrimer is detected in the current
        orientation) the read will either be written to the forward (read 1) or
        reverse (read 2) reads for the case of paired files, or the read will be
        reverse complemented in the case of stitched reads.
    disable_header_match: if True, suppresses checks between fastq headers.
    """

    # Turn off extra file creation for single read.
    if input_type == "barcode_single_end" and attempt_read_orientation:
        attempt_read_orientation = False
    if attempt_read_orientation:
        header, mapping_data, run_description, errors, warnings =\
            process_id_map(map_fp)
        forward_primers, reverse_primers = get_primers(header, mapping_data)
        output_bc_not_oriented = open(join(output_dir,
                                           "barcodes_not_oriented.fastq.incomplete"), "w")
        fastq1_out_not_oriented = open(join(output_dir,
                                            "reads1_not_oriented.fastq.incomplete"), "w")
        fastq2_out_not_oriented = open(join(output_dir,
                                            "reads2_not_oriented.fastq.incomplete"), "w")
    else:
        forward_primers = None
        reverse_primers = None
        output_bc_not_oriented = None
        fastq1_out_not_oriented = None
        fastq2_out_not_oriented = None

    output_bc_fastq = open(join(output_dir, "barcodes.fastq.incomplete"), "w")
    if input_type in ["barcode_single_end", "barcode_paired_stitched"]:
        output_fastq1 = open(join(output_dir, "reads.fastq.incomplete"), "w")
        output_fastq2 = None
        final_fastq1_name = join(output_dir, "reads.fastq")
    elif input_type in ["barcode_paired_end"]:
        output_fastq1 = open(join(output_dir, "reads1.fastq.incomplete"), "w")
        output_fastq2 = open(join(output_dir, "reads2.fastq.incomplete"), "w")
        final_fastq1_name = join(output_dir, "reads1.fastq")
    else:
        output_fastq1 = None
        output_fastq2 = None

    if not fastq2:
        fastq2 = cycle(["@", "AAAAAAAAAAAA", "+", "AAAAAAAAAAAA"])
        not_paired = True
    else:
        not_paired = False

    check_header_match_f = get_casava_version(fastq1)

    header_index = 0

    for read1_data, read2_data in izip(
            parse_fastq(fastq1, strict=False, enforce_qual_range=False),
            parse_fastq(fastq2, strict=False, enforce_qual_range=False)):
        if not disable_header_match:
            if not check_header_match_f(read1_data[header_index],
                                        read2_data[header_index]):
                raise FastqParseError("Headers of read1 and read2 do not match. Can't continue. "
                                      "Confirm that the fastq sequences that you are "
                                      "passing match one another. --disable_header_match can be "
                                      "used to suppress header checks.")

        if input_type == "barcode_single_end":
            process_barcode_single_end_data(read1_data, output_bc_fastq,
                                            output_fastq1, bc1_len, rev_comp_bc1)

        elif input_type == "barcode_paired_end":
            process_barcode_paired_end_data(read1_data, read2_data,
                                            output_bc_fastq, output_fastq1, output_fastq2, bc1_len, bc2_len,
                                            rev_comp_bc1, rev_comp_bc2, attempt_read_orientation,
                                            forward_primers, reverse_primers, output_bc_not_oriented,
                                            fastq1_out_not_oriented, fastq2_out_not_oriented)

        elif input_type == "barcode_paired_stitched":
            process_barcode_paired_stitched(read1_data,
                                            output_bc_fastq, output_fastq1, bc1_len, bc2_len,
                                            rev_comp_bc1, rev_comp_bc2, attempt_read_orientation,
                                            forward_primers, reverse_primers, output_bc_not_oriented,
                                            fastq1_out_not_oriented, switch_bc_order)

        elif input_type == "barcode_in_label":
            if not_paired:
                curr_read2_data = False
            else:
                curr_read2_data = read2_data
            process_barcode_in_label(read1_data, curr_read2_data,
                                     output_bc_fastq, bc1_len, bc2_len,
                                     rev_comp_bc1, rev_comp_bc2, char_delineator)

    output_bc_fastq.close()
    rename(output_bc_fastq.name, join(output_dir, "barcodes.fastq"))
    if output_fastq1:
        output_fastq1.close()
        rename(output_fastq1.name, final_fastq1_name)
    if output_fastq2:
        output_fastq2.close()
        rename(output_fastq2.name, join(output_dir, "reads2.fastq"))
    if output_bc_not_oriented:
        rename(output_bc_not_oriented.name,
               join(output_dir, "barcodes_not_oriented.fastq"))
    if fastq1_out_not_oriented:
        rename(fastq1_out_not_oriented.name,
               join(output_dir, "reads1_not_oriented.fastq"))
    if fastq2_out_not_oriented:
        rename(fastq2_out_not_oriented.name,
               join(output_dir, "reads2_not_oriented.fastq"))


def get_casava_version(fastq1):
    """ determine the version of casava that was used to generate the fastq

    fastq1: open fastq file or list of lines
    """
    # to determine how to compare header lines and decode ascii phred scores
    # Modified code from split_libraries_fastq.py for this check
    if isinstance(fastq1, list):
        fastq_read_f_line1 = fastq1[0]
    else:
        fastq_read_f_line1 = fastq1.readline()
        fastq1.seek(0)
    post_casava_v180 = is_casava_v180_or_later(fastq_read_f_line1)
    if post_casava_v180:
        check_header_match_f = check_header_match_180_or_later
    else:
        check_header_match_f = check_header_match_pre180

    return check_header_match_f


def process_barcode_single_end_data(read1_data,
                                    output_bc_fastq,
                                    output_fastq1,
                                    bc1_len=6,
                                    rev_comp_bc1=False):
    """ Processes, writes single-end barcode data, parsed sequence

    read1_data: list of header, read, quality scores
    output_bc_fastq: open output fastq filepath
    output_fastq1: open output fastq reads filepath
    bc1_len: length of barcode to remove from beginning of data
    rev_comp_bc1: reverse complement barcode before writing.
    """

    header_index = 0
    sequence_index = 1
    quality_index = 2

    bc_read = read1_data[sequence_index][:bc1_len]
    bc_qual = read1_data[quality_index][:bc1_len]
    if rev_comp_bc1:
        bc_read = str(DNA(bc_read).rc())
        bc_qual = bc_qual[::-1]

    bc_lines = format_fastq_record(read1_data[header_index], bc_read, bc_qual)
    output_bc_fastq.write(bc_lines)
    seq_lines = format_fastq_record(read1_data[header_index],
                                    read1_data[sequence_index][bc1_len:],
                                    read1_data[quality_index][bc1_len:])
    output_fastq1.write(seq_lines)

    return


def process_barcode_paired_end_data(read1_data,
                                    read2_data,
                                    output_bc_fastq,
                                    output_fastq1,
                                    output_fastq2,
                                    bc1_len=6,
                                    bc2_len=6,
                                    rev_comp_bc1=False,
                                    rev_comp_bc2=False,
                                    attempt_read_orientation=False,
                                    forward_primers=None,
                                    reverse_primers=None,
                                    output_bc_not_oriented=None,
                                    fastq1_out_not_oriented=None,
                                    fastq2_out_not_oriented=None):
    """ Processes, writes paired-end barcode data, parsed sequences

    read1_data: list of header, read, quality scores
    read2_data: list of header, read, quality scores
    output_bc_fastq: open output fastq filepath
    output_fastq1: open output fastq reads 1 filepath
    output_fastq2: open output fastq reads 2 filepath
    bc1_len: length of barcode to remove from beginning of read1 data
    bc2_len: length of barcode to remove from beginning of read2 data
    rev_comp_bc1: reverse complement barcode 1 before writing.
    rev_comp_bc2: reverse complement barcode 2 before writing.
    attempt_read_orientation: If True, will attempt to orient the reads
        according to the forward primers in the mapping file. If primer is
        detected in current orientation, leave the read as is, but if reverse
        complement is detected (or ReversePrimer is detected in the current
        orientation) the read will either be written to the forward (read 1) or
        reverse (read 2) reads for the case of paired files, or the read will be
        reverse complemented in the case of stitched reads.
    forward_primers: list of regular expression generators, forward primers
    reverse_primers: list of regular expression generators, reverse primers
    output_bc_not_oriented: Barcode output from reads that are not oriented
    fastq1_out_not_oriented: Open filepath to write reads 1 where primers
        can't be found when attempt_read_orientation is True.
    fastq2_out_not_oriented: Open filepath to write reads 2 where primers
        can't be found when attempt_read_orientation is True.
    """

    header_index = 0
    sequence_index = 1
    quality_index = 2

    found_primer_match = False
    # Break from orientation search as soon as a match is found
    if attempt_read_orientation:
        # First check forward primers
        for curr_primer in forward_primers:
            if curr_primer.search(read1_data[sequence_index]):
                read1 = read1_data
                read2 = read2_data
                found_primer_match = True
                break
            if curr_primer.search(read2_data[sequence_index]):
                read1 = read2_data
                read2 = read1_data
                found_primer_match = True
                break
        # Check reverse primers if forward primers not found
        if not found_primer_match:
            for curr_primer in reverse_primers:
                if curr_primer.search(read1_data[sequence_index]):
                    read1 = read2_data
                    read2 = read1_data
                    found_primer_match = True
                    break
                if curr_primer.search(read2_data[sequence_index]):
                    read1 = read1_data
                    read2 = read2_data
                    found_primer_match = True
                    break
    else:
        read1 = read1_data
        read2 = read2_data

    if not found_primer_match and attempt_read_orientation:
        read1 = read1_data
        read2 = read2_data
        output_bc = output_bc_not_oriented
        output_read1 = fastq1_out_not_oriented
        output_read2 = fastq2_out_not_oriented
    else:
        output_bc = output_bc_fastq
        output_read1 = output_fastq1
        output_read2 = output_fastq2

    bc_read1 = read1[sequence_index][0:bc1_len]
    bc_read2 = read2[sequence_index][0:bc2_len]
    bc_qual1 = read1[quality_index][0:bc1_len]
    bc_qual2 = read2[quality_index][0:bc2_len]
    if rev_comp_bc1:
        bc_read1 = str(DNA(bc_read1).rc())
        bc_qual1 = bc_qual1[::-1]
    if rev_comp_bc2:
        bc_read2 = str(DNA(bc_read2).rc())
        bc_qual2 = bc_qual2[::-1]

    bc_lines = format_fastq_record(read1[header_index],
                                   bc_read1 + bc_read2,
                                   np.hstack([bc_qual1, bc_qual2]))
    output_bc.write(bc_lines)
    seq1_lines = format_fastq_record(read1[header_index],
                                     read1[sequence_index][bc1_len:], read1[quality_index][bc1_len:])
    output_read1.write(seq1_lines)
    seq2_lines = format_fastq_record(read2[header_index],
                                     read2[sequence_index][bc2_len:], read2[quality_index][bc2_len:])
    output_read2.write(seq2_lines)

    return


def process_barcode_paired_stitched(read_data,
                                    output_bc_fastq,
                                    output_fastq,
                                    bc1_len=6,
                                    bc2_len=6,
                                    rev_comp_bc1=False,
                                    rev_comp_bc2=False,
                                    attempt_read_orientation=False,
                                    forward_primers=None,
                                    reverse_primers=None,
                                    output_bc_not_oriented=None,
                                    fastq_out_not_oriented=None,
                                    switch_bc_order=False):
    """ Processes stitched barcoded reads, writes barcode, parsed stitched read

    read_data: list of header, read, quality scores
    output_bc_fastq: open output fastq filepath
    output_fastq: open output fastq reads filepath
    bc1_len: length of barcode to remove from beginning of read1 stitched data
    bc2_len: length of barcode to remove from end of read2 stitched data
    rev_comp_bc1: reverse complement barcode 1 before writing.
    rev_comp_bc2: reverse complement barcode 2 before writing.
    attempt_read_orientation: If True, will attempt to orient the reads
        according to the forward primers in the mapping file. If primer is
        detected in current orientation, leave the read as is, but if reverse
        complement is detected (or ReversePrimer is detected in the current
        orientation) the read will either be written to the forward (read 1) or
        reverse (read 2) reads for the case of paired files, or the read will be
        reverse complemented in the case of stitched reads.
    forward_primers: list of regular expression generators, forward primers
    reverse_primers: list of regular expression generators, reverse primers
    output_bc_not_oriented: Barcode output from reads that are not oriented
    fastq_out_not_oriented: Open filepath to write reads where primers
        can't be found when attempt_read_orientation is True.
    switch_bc_order: Normally, barcode 1 will be written first, followed by
        barcode 2 in a combined output fastq file. If True, the order will be
        reversed. Only applies to stitched reads processing, as other barcode
        orders are dictated by the the parameter chosen for the fastq files.
    """

    header_index = 0
    sequence_index = 1
    quality_index = 2

    read_seq = read_data[sequence_index]
    read_qual = read_data[quality_index]

    found_primer_match = False
    # Break from orientation search as soon as a match is found
    if attempt_read_orientation:
        for curr_primer in forward_primers:
            if curr_primer.search(read_data[sequence_index]):
                found_primer_match = True
                break
        if not found_primer_match:
            for curr_primer in reverse_primers:
                if curr_primer.search(read_data[sequence_index]):
                    read_seq = str(DNA(read_seq).rc())
                    read_qual = read_qual[::-1]
                    found_primer_match = True
                    break

    if not found_primer_match and attempt_read_orientation:
        output_bc = output_bc_not_oriented
        output_read = fastq_out_not_oriented
    else:
        output_bc = output_bc_fastq
        output_read = output_fastq

    bc_read1 = read_seq[0:bc1_len]
    bc_read2 = read_seq[-bc2_len:]
    bc_qual1 = read_qual[0:bc1_len]
    bc_qual2 = read_qual[-bc2_len:]

    if rev_comp_bc1:
        bc_read1 = str(DNA(bc_read1).rc())
        bc_qual1 = bc_qual1[::-1]
    if rev_comp_bc2:
        bc_read2 = str(DNA(bc_read2).rc())
        bc_qual2 = bc_qual2[::-1]

    if switch_bc_order:
        bc_read1, bc_read2 = bc_read2, bc_read1
        bc_qual1, bc_qual2 = bc_qual2, bc_qual1

    bc_lines = format_fastq_record(read_data[header_index],
                                   bc_read1 + bc_read2,
                                   np.hstack([bc_qual1, bc_qual2]))
    output_bc.write(bc_lines)
    seq_lines = format_fastq_record(read_data[header_index],
                                    read_seq[bc1_len:-bc2_len], read_qual[bc1_len:-bc2_len])
    output_read.write(seq_lines)

    return


def process_barcode_in_label(read1_data,
                             read2_data,
                             output_bc_fastq,
                             bc1_len=6,
                             bc2_len=6,
                             rev_comp_bc1=False,
                             rev_comp_bc2=False,
                             char_delineator=":"):
    """ Reads data from one or two fastq labels, writes output barcodes file.

    read1_data: list of header, read, quality scores
    read2_data: list of header, read, quality scores, False if no read 2.
    output_bc_fastq: open output fastq filepath
    bc1_len: length of barcode to remove from beginning of read1 data
    bc2_len: length of barcode to remove from beginning of read2 data
    rev_comp_bc1: reverse complement barcode 1 before writing.
    rev_comp_bc2: reverse complement barcode 2 before writing.
    char_delineator: Specify character that immediately precedes the barcode
        for input_type of barcode_in_label.
    """
    header_index = 0

    # Check for char_delineator in sequence
    try:
        bc1_read = read1_data[header_index].split(
            char_delineator)[-1][0:bc1_len]
    # If there is an index error, it means the char_delineator wasn't found
    except IndexError:
        raise IndexError("Found sequence lacking character delineator. "
                         "Sequence header %s, character delineator %s" %
                         (read1_data[header_index], char_delineator))

    # Create fake quality scores, using 6 here to match the existing qual fake
    # qual scores that were all F.
    bc1_qual = np.ones(len(bc1_read), dtype=np.int8) * 6
    if rev_comp_bc1:
        bc1_read = str(DNA(bc1_read).rc())

    if read2_data:
        bc2_read =\
            read2_data[header_index].strip().split(
                char_delineator)[-1][0:bc2_len]
        bc2_qual = np.ones(len(bc2_read), dtype=np.int8) * 6
        if rev_comp_bc2:
            bc2_read = str(DNA(bc2_read).rc())
    else:
        bc2_read = ""
        bc2_qual = np.array([], dtype=np.int8)

    if not bc1_read and not bc2_read:
        raise ValueError("Came up with empty barcode sequence, please check "
                         "character delineator with -s, and fastq label "
                         "%s" % read1_data[header_index])

    bc_lines = format_fastq_record(read1_data[header_index],
                                   bc1_read + bc2_read,
                                   np.hstack([bc1_qual, bc2_qual]))

    output_bc_fastq.write(bc_lines)

    return


def get_primers(header,
                mapping_data):
    """ Returns lists of forward/reverse primer regular expression generators

    header:  list of strings of header data.
    mapping_data:  list of lists of mapping data

    Will raise error if either the LinkerPrimerSequence or ReversePrimer fields
        are not present
    """

    if "LinkerPrimerSequence" in header:
        primer_ix = header.index("LinkerPrimerSequence")
    else:
        raise IndexError(
            ("Mapping file is missing LinkerPrimerSequence field."))
    if "ReversePrimer" in header:
        rev_primer_ix = header.index("ReversePrimer")
    else:
        raise IndexError(("Mapping file is missing ReversePrimer field."))

    iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
             'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
             'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}

    raw_forward_primers = set([])
    raw_forward_rc_primers = set([])
    raw_reverse_primers = set([])
    raw_reverse_rc_primers = set([])

    for line in mapping_data:
        # Split on commas to handle pool of primers
        raw_forward_primers.update([upper(primer).strip() for
                                    primer in line[primer_ix].split(',')])
        raw_forward_rc_primers.update([str(DNA(primer).rc()) for
                                       primer in raw_forward_primers])
        raw_reverse_primers.update([upper(primer).strip() for
                                    primer in line[rev_primer_ix].split(',')])
        raw_reverse_rc_primers.update([str(DNA(primer).rc()) for
                                       primer in raw_reverse_primers])

    if not raw_forward_primers:
        raise ValueError(("No forward primers detected in mapping file."))
    if not raw_reverse_primers:
        raise ValueError(("No reverse primers detected in mapping file."))

    # Finding the forward primers, or rc of reverse primers indicates forward
    # read. Finding the reverse primer, or rc of the forward primers, indicates
    # the reverse read, so these sets are merged.
    raw_forward_primers.update(raw_reverse_rc_primers)
    raw_reverse_primers.update(raw_forward_rc_primers)

    forward_primers = []
    reverse_primers = []
    for curr_primer in raw_forward_primers:
        forward_primers.append(compile(''.join([iupac[symbol] for
                                                symbol in curr_primer])))
    for curr_primer in raw_reverse_primers:
        reverse_primers.append(compile(''.join([iupac[symbol] for
                                                symbol in curr_primer])))

    return forward_primers, reverse_primers
