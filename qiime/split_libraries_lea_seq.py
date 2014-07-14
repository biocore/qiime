#!/usr/bin/env python
from __future__ import division

__author__ = "Charudatta Navare"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Charudatta Navare", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Charudatta Navare"
__email__ = "charudatta.navare@gmail.com"

from collections import defaultdict
from qiime.golay import get_invalid_golay_barcodes
from qiime.parse import parse_mapping_file_to_dict
from qiime.split_libraries import check_map, expand_degeneracies
from qiime.split_libraries_fastq import correct_barcode
from skbio.parse.sequences import parse_fasta, parse_fastq
from qiime.util import get_qiime_temp_dir, qiime_system_call
from brokit.uclust import get_clusters_from_fasta_filepath
from skbio.app.util import ApplicationError
from tempfile import mkstemp
from itertools import izip
from os import close, path, mkdir, rmdir
from skbio.util.misc import remove_files
from skbio.core.sequence import DNASequence
import re


class PairedEndParseError(Exception):
    pass


class PrimerMismatchError(Exception):
    pass


class LowConsensusScoreError(Exception):
    pass


class SeqLengthMismatchError(Exception):
    pass


class InvalidGolayBarcodeError(Exception):
    pass


class BarcodeLenMismatchError(Exception):
    pass


def extract_primer(seq, possible_primers, min_idx=None, max_idx=None):
    """
    Extracts primers from sequence, given possible primers
    returns before_primer, primer, sequence without primer

    """

    primer_idx = None   # index of primer in sequence
    primer = None      # primer in sequence

    for possible_primer in possible_primers:
        if possible_primer in seq:
            primer_idx = seq.index(possible_primer)
            primer = possible_primer

            if(min_idx is not None and primer_idx < min_idx) or \
              (max_idx is not None and primer_idx > max_idx):
                primer_idx = None
                primer = None
                continue
            else:
                break

    if primer_idx is None:
        raise PrimerMismatchError(
            "Sequence does not contain any primers from the mapping file, "
            "Please verify your mapping file.")

    before_primer = seq[:primer_idx]
    after_primer = seq.replace(before_primer + primer, '', 1)
    return before_primer, primer, after_primer


def get_LEA_seq_consensus_seqs(fwd_read_f, rev_read_f,
                               map_f, output_dir,
                               barcode_type, barcode_len,
                               barcode_correction_fn, max_barcode_errors,
                               min_consensus, max_cluster_ratio,
                               min_difference_in_bcs,
                               fwd_length, rev_length,
                               min_reads_per_random_bc,
                               min_difference_in_clusters,
                               barcode_column,
                               reverse_primer_column):
    """
    Reads mapping file, input file, and other command line arguments
    fills dictionary called consensus_seq_lookup which will contain:
    sample ID -> random barcode -> consensus_seq
    """

    (bc_to_sid,
     bc_to_fwd_primers,
     bc_to_rev_primers) = process_mapping_file(map_f,
                                               barcode_len,
                                               barcode_type,
                                               barcode_column,
                                               reverse_primer_column)

    (random_bc_lookup,
     random_bc_reads,
     random_bcs,
     barcode_errors_exceed_max_count,
     barcode_not_in_map_count,
     primer_mismatch_count,
     seq_too_short_count,
     input_seqs_count,
     total_seqs_kept) = read_fwd_rev_read(fwd_read_f,
                                          rev_read_f,
                                          bc_to_sid,
                                          barcode_len,
                                          barcode_correction_fn,
                                          bc_to_fwd_primers,
                                          bc_to_rev_primers,
                                          max_barcode_errors,
                                          fwd_length,
                                          rev_length)

    consensus_seq_lookup = get_consensus_seqs_lookup(random_bc_lookup,
                                                     random_bc_reads,
                                                     random_bcs,
                                                     min_difference_in_bcs,
                                                     min_reads_per_random_bc,
                                                     output_dir,
                                                     min_difference_in_clusters,
                                                     max_cluster_ratio)

    log_out = format_lea_seq_log(input_seqs_count,
                                 barcode_errors_exceed_max_count,
                                 barcode_not_in_map_count,
                                 primer_mismatch_count,
                                 seq_too_short_count,
                                 total_seqs_kept)

    return consensus_seq_lookup, log_out


def get_cluster_ratio(fasta_seqs, min_difference_in_clusters):
    """
    Uses uclust to calculate cluster ratio
    cluster_ratio=num_of_seq_in_cluster_with_max_seq/num_of_seq_in cluster_with_second_higest_seq
    """
    cluster_percent_id = min_difference_in_clusters
    temp_dir = get_qiime_temp_dir()
    fd_uc, uclust_tempfile_name = mkstemp(dir=temp_dir, suffix='.uc')
    close(fd_uc)
    fd_fas, fasta_tempfile_name = mkstemp(dir=temp_dir, suffix='.uc')
    close(fd_fas)
    fasta_tempfile = open(fasta_tempfile_name, 'w')
    fasta_tempfile.write(fasta_seqs)
    fasta_tempfile.close()
    count_lookup = {}
    count = 0
    command = "uclust --usersort --input " + fasta_tempfile_name +\
              " --uc " + uclust_tempfile_name + " --id 0.98"
    # In the function, I am calling uclust a large number of times.
    # Initially I was using brokit.get_clusters_from_fasta_filepath
    # but due to issue (biocore/brokit#31), I have temporarily
    # reverted to qiime_system_call.

    qiime_system_call(command)
    uclust_tempfile = open(uclust_tempfile_name, 'r')
    for line in uclust_tempfile:
        if re.search(r'^C', line):
            pieces = line.split('\t')
            try:
                count_lookup[pieces[1]] += pieces[2]
            except KeyError:
                count_lookup[pieces[1]] = pieces[2]
            except IndexError:
                pass
            count += 1
    uclust_tempfile.close()
    files_to_be_removed = list()
    files_to_be_removed.append(uclust_tempfile_name)
    remove_files(files_to_be_removed)

    sorted_counts_in_clusters = sorted(
        count_lookup.iteritems(),
        key=lambda x: x[1])
    try:
        max_cluster_count = \
            float(str(sorted_counts_in_clusters[0][1]))
        second_cluster_count = \
            float(str(sorted_counts_in_clusters[1][1]))
        return max_cluster_count / second_cluster_count
    except IndexError:
        return 1


def get_consensus(fasta_seqs, min_consensus):
    """
    Returns consensus sequence from a set of sequences
    input: fasta file, min_consensus
    fasta_file should be in the following format:
    >random_bc|number
    seq
    >random_bc|number
    seq
    ....

    number = number of times this seq has appeared with this random_barcode
    """
    seqs = list()
    counts = list()

    temp_dir = get_qiime_temp_dir()
    fd_fas, fasta_tempfile_name = mkstemp(dir=temp_dir, suffix='.fas')
    close(fd_fas)

    fasta_tempfile = open(fasta_tempfile_name, 'w')
    fasta_tempfile.write(fasta_seqs)
    fasta_tempfile.close()

    fasta_tempfile = open(fasta_tempfile_name, 'r')
    for label, seq in parse_fasta(fasta_tempfile):
        RE_output = re.search(r'\w+\|(\d+)', label)
        counts.append(int(RE_output.group(1)))
        seqs.append(seq)
    fasta_tempfile.close()

    files_to_be_removed = list()
    files_to_be_removed.append(fasta_tempfile_name)
    remove_files(files_to_be_removed)

    length = len(seqs[0])
    number_of_seqs = len(seqs)

    for seq_index in range(number_of_seqs):
        if len(seqs[seq_index]) != length:
            raise SeqLengthMismatchError()

    freq_this_pos_this_base = dict()
    count_of_seq_with_max_count = dict()

    for x in range(length):
        freq_this_pos_this_base[x] = dict()
        count_of_seq_with_max_count[x] = dict()
    for x in range(length):
        for y in DNASequence.iupac_characters():
            freq_this_pos_this_base[x][y] = 0
            count_of_seq_with_max_count[x][y] = 0

    for base_index in range(length):
        for this_seq_count, seq in enumerate(seqs):
            freq_this_pos_this_base[base_index][
                seq[base_index]] += counts[this_seq_count]
            if counts[this_seq_count] > count_of_seq_with_max_count[
                    base_index][seq[base_index]]:
                count_of_seq_with_max_count[base_index][
                    seq[base_index]] = counts[this_seq_count]

    consensus = list()
    for index in range(length):
        sorted_bases = sorted(
            freq_this_pos_this_base[index].iteritems(),
            key=lambda x: x[1])
        max_base, max_freq = sorted_bases[-1]

        for (counter, (b, n)) in enumerate(sorted_bases):
            if max_freq == n:
                try:
                    if count_of_seq_with_max_count[counter][
                            b] > count_of_seq_with_max_count[counter][max_base]:
                        max_base = b
                except KeyError:
                    pass

        score = 10.0 * max_freq / number_of_seqs
        if score < min_consensus:
            raise LowConsensusScoreError()
        consensus.append(max_base)

    consensus_seq = ''.join(map(str, consensus))
    return consensus_seq


def select_unique_rand_bcs(rand_bcs, min_difference_in_bcs):
    """
    Attempts to select true barcodes from set of barcodes
    i.e. removes barcodes that might be artifacts
    due to sequencing errors.
    Uses uclust to remove barcodes that are similar thatn
    threshold.
    returns: a set containing random unique random barcodes.
    """
    unique_threshold = min_difference_in_bcs
    temp_dir = get_qiime_temp_dir()
    fasta_fd, fasta_tempfile_name = mkstemp(
        dir=temp_dir, prefix='tmp', suffix='.fas')
    rand_bcs = set(rand_bcs)
    fasta_tempfile = open(fasta_tempfile_name, "w")
    p_line = ""
    for rand_bc in rand_bcs:
        p_line = p_line + ">" + rand_bc + "\n" + rand_bc + "\n"
    fasta_tempfile.write(p_line)
    fasta_tempfile.close()

    _, _, unique_rand_bcs = get_clusters_from_fasta_filepath(
        fasta_tempfile_name,
        original_fasta_path=None,
        percent_ID=unique_threshold,
        save_uc_files=False,
        output_dir=temp_dir)

    unique_rand_bcs = set(unique_rand_bcs)

    files_to_be_removed = list()
    files_to_be_removed.append(fasta_tempfile_name)
    remove_files(files_to_be_removed)
    return unique_rand_bcs


def format_lea_seq_log(input_seqs_count,
                       barcode_errors_exceed_max_count,
                       barcode_not_in_map_count,
                       primer_mismatch_count,
                       seq_too_short_count,
                       total_seqs_kept):
    """ Format the split libraries LEA-Seq log """
    log_out = ["Quality filter results"]
    log_out.append(
        "Total number of input sequences: %d" %
        input_seqs_count)
    log_out.append(
        "Barcode not in mapping file: %d" %
        barcode_not_in_map_count)
    log_out.append(
        "Sequence shorter than threshold: %d" %
        seq_too_short_count)
    log_out.append(
        "Barcode errors exceeds limit: %d" %
        barcode_errors_exceed_max_count)
    log_out.append(
        "Primer mismatch count: %d" %
        primer_mismatch_count)
    log_out.append("")

    log_out.append('\nTotal number seqs written: %d' % total_seqs_kept)
    return '\n'.join(log_out)


def check_barcodes(bc_to_sid, barcode_len, barcode_type):
    """
    Make sure that barcodes (which are guaranteed to be of
    the same length at this point) are the correct length
    that the user specified.
    """
    barcode_len_in_map = len(bc_to_sid.keys()[0])
    if barcode_len_in_map != barcode_len:
        raise BarcodeLenMismatchError("Barcodes in mapping file are of length %d, but "
                                      "expected barcodes of length %d." %
                                      (barcode_len_in_map, barcode_len))

    if barcode_type == 'golay_12':
        invalid_golay_barcodes = get_invalid_golay_barcodes(bc_to_sid.keys())

        if invalid_golay_barcodes:
            raise InvalidGolayBarcodeError(
                "Some or all barcodes in the mapping file are "
                "not valid golay codes. Do they need to be "
                "reverse complemented? If these are not golay "
                "barcodes pass --barcode_type 12 to disable "
                "barcode error correction, or pass "
                "--barcode_type  # if the barcodes are not 12 "
                "base pairs, where   #  is the size of the "
                "barcodes.\n\nInvalid barcodes: %s" %
                ' '.join(invalid_golay_barcodes))


def get_consensus_seqs_lookup(random_bc_lookup,
                              random_bc_reads,
                              random_bcs,
                              min_difference_in_bcs,
                              min_reads_per_random_bc,
                              output_dir,
                              min_difference_in_clusters,
                              max_cluster_ratio):
    """
    Generate LEA-seq consensus sequence
    returns defaultdict called consensus_seq_lookup
    which will contain:
    sample ID -> random barcode -> consensus_seq
    """

    consensus_seq_lookup = defaultdict(lambda:
                                       defaultdict(str))
    random_bc_keep = {}

    for sample_id in random_bc_lookup:
        random_bc_keep[sample_id] = select_unique_rand_bcs(
            random_bcs[sample_id],
            min_difference_in_bcs)
        for random_bc in random_bc_lookup[sample_id]:
            if random_bc in random_bc_keep[sample_id] and random_bc_reads[
                    sample_id][random_bc] >= min_reads_per_random_bc:
                fwd_fd, fwd_fasta_tempfile_name = mkstemp(
                    dir=output_dir, prefix='fwd', suffix='.fas')
                rev_fd, rev_fasta_tempfile_name = mkstemp(
                    dir=output_dir, prefix='rev', suffix='.fas')
                close(fwd_fd)
                close(rev_fd)
                fwd_fasta_tempfile = open(fwd_fasta_tempfile_name, 'w')
                rev_fasta_tempfile = open(rev_fasta_tempfile_name, 'w')
                max_freq = 0
                for seq_index, fwd_rev in enumerate(random_bc_lookup[sample_id][random_bc]):
                    fwd_seq, rev_seq = fwd_rev
                    fwd_line = ">" + str(seq_index) + random_bc + "|" + str(
                        random_bc_lookup[sample_id][random_bc][fwd_rev]) +\
                        "\n" + fwd_seq + "\n"
                    rev_line = ">" + str(seq_index) + random_bc + "|" + str(
                        random_bc_lookup[sample_id][random_bc][fwd_rev]) +\
                        "\n" + rev_seq + "\n"
                    fwd_fasta_tempfile.write(fwd_line)
                    rev_fasta_tempfile.write(rev_line)
                    if random_bc_lookup[sample_id][
                            random_bc][fwd_rev] > max_freq:
                        max_freq = random_bc_lookup[
                            sample_id][random_bc][fwd_rev]
                        majority_seq = fwd_seq + "^" + rev_seq
                fwd_fasta_tempfile.close()
                rev_fasta_tempfile.close()
                fwd_cluster_ratio = get_cluster_ratio(
                    fwd_fasta_tempfile_name,
                    min_difference_in_clusters)
                rev_cluster_ratio = get_cluster_ratio(
                    rev_fasta_tempfile_name,
                    min_difference_in_clusters)
                if fwd_cluster_ratio == 0 or rev_cluster_ratio == 0:
                    consensus_seq = "No consensus"
                elif fwd_cluster_ratio < max_cluster_ratio and rev_cluster_ratio < max_cluster_ratio:
                    consensus_seq = majority_seq
                else:
                    fwd_fasta_tempfile = open(fwd_fasta_tempfile_name, 'r')
                    rev_fasta_tempfile = open(rev_fasta_tempfile_name, 'r')
                    fwd_consensus = get_consensus(
                        fwd_fasta_tempfile,
                        min_consensus)
                    rev_consensus = get_consensus(
                        rev_fasta_tempfile,
                        min_consensus)
                    fwd_fasta_tempfile.close()
                    rev_fasta_tempfile.close()
                    consensus_seq = fwd_consensus + "^" + rev_consensus

                consensus_seq_lookup[sample_id][random_bc] = consensus_seq
                files_to_be_removed = list()
                files_to_be_removed.append(fwd_fasta_tempfile_name)
                files_to_be_removed.append(rev_fasta_tempfile_name)
                remove_files(files_to_be_removed)
                return consensus_seq_lookup


def read_fwd_rev_read(fwd_read_f,
                      rev_read_f,
                      bc_to_sid,
                      barcode_len,
                      barcode_correction_fn,
                      bc_to_fwd_primers,
                      bc_to_rev_primers,
                      max_barcode_errors,
                      fwd_length,
                      rev_length):

    random_bc_lookup = defaultdict(lambda:
                                   defaultdict(lambda:
                                               defaultdict(int)))

    random_bc_reads = defaultdict(lambda:
                                  defaultdict(int))

    random_bcs = {}

    # Counts for Quality Control:
    input_seqs_count = 0
    total_seqs_kept_count = 0
    barcode_errors_exceed_max_count = 0
    barcode_not_in_map_count = 0
    primer_mismatch_count = 0
    seq_too_short_count = 0
    input_seqs_count = 0
    total_seqs_kept = 0

    header_idx = 0
    seq_idx = 1
    qual_idx = 2

    for fwd_read, rev_read in izip(
            parse_fastq(fwd_read_f, strict=False),
            parse_fastq(rev_read_f, strict=False)):
            # Confirm match between read headers.

        input_seqs_count += 1

        if fwd_read[header_idx] != rev_read[header_idx]:
            raise PairedEndParseError(
                "Headers of forward and reverse reads "
                "do not match. Confirm that the forward "
                "and reverse read fastq files that you "
                "provided have headers that match one "
                "another.")
        else:
            header = fwd_read[header_idx]

        fwd_seq = fwd_read[seq_idx]
        rev_seq = rev_read[seq_idx]

        #  Grab the barcode sequence. It is always at the very end of the
        #  forward read. Strip the barcode from the sequence.
        barcode = fwd_seq[-barcode_len:]
        fwd_seq = fwd_seq[:-barcode_len]

        #  Correct the barcode(if applicable) and map to sample ID.
        num_barcode_errors, corrected_barcode, _, sample_id =\
            correct_barcode(barcode, bc_to_sid, barcode_correction_fn)

        #  Skip barcodes with too many errors.
        if num_barcode_errors > max_barcode_errors:
            barcode_errors_exceed_max_count += 1
            continue

        if sample_id is None:
            barcode_not_in_map_count += 1
            continue

        # Extract the random barcode and primer from the forward read.
        possible_primers = bc_to_fwd_primers[corrected_barcode].keys()

        try:
            random_bc, _, clean_fwd_seq = extract_primer(fwd_seq,
                                                         possible_primers,
                                                         min_idx=5,
                                                         max_idx=20)

            random_bcs[sample_id].append(random_bc)
        except PrimerMismatchError:
            primer_mismatch_count += 1
            continue
        except KeyError:
            random_bcs[sample_id] = list()
            random_bcs[sample_id].append(random_bc)

        possible_primers = bc_to_rev_primers[barcode]

        try:
            phase_seq, _, clean_rev_seq = extract_primer(rev_seq,
                                                         possible_primers)
        except PrimerMismatchError:
            primer_mismatch_count += 1
            continue

        if len(clean_fwd_seq) < fwd_length:
            seq_too_short_count += 1
            continue

        clean_fwd_seq = clean_fwd_seq[:fwd_length]
        clean_rev_seq = clean_rev_seq[:rev_length]

        total_seqs_kept += 1
        random_bc_reads[sample_id][random_bc] += 1
        random_bc_lookup[sample_id][random_bc][
            (clean_fwd_seq, clean_rev_seq)] += 1

    return (random_bc_lookup,
            random_bc_reads,
            random_bcs,
            barcode_errors_exceed_max_count,
            barcode_not_in_map_count,
            primer_mismatch_count,
            seq_too_short_count,
            input_seqs_count,
            total_seqs_kept)


def process_mapping_file(map_f,
                         barcode_len,
                         barcode_type,
                         BARCODE_COLUMN,
                         REVERSE_PRIMER_COLUMN):
    """Ensures that sample IDs and barcodes are unique, that barcodes are
    all the same length, and that primers are present. Ensures barcodes
    and primers only contain valid characters."""

    _, _, bc_to_sid, _, _, bc_to_fwd_primers, _ = check_map(map_f, False)
    map_f.seek(0)

    metadata_map = parse_mapping_file_to_dict(map_f)[0]
    bc_to_rev_primers = {}
    for sid, md in metadata_map.items():
        if REVERSE_PRIMER_COLUMN in md:
            bc_to_rev_primers[
                md[BARCODE_COLUMN]] = expand_degeneracies(
                md[REVERSE_PRIMER_COLUMN].upper().split(','))
        else:
            raise Exception(
                "The %s column does not exist in the "
                "mapping file. %s is required." %
                (REVERSE_PRIMER_COLUMN,
                 REVERSE_PRIMER_COLUMN))

    check_barcodes(bc_to_sid, barcode_len, barcode_type)

    return (bc_to_sid,
            bc_to_fwd_primers,
            bc_to_rev_primers)
