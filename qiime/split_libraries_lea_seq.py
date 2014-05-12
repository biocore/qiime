from collections import defaultdict
from optparse import OptionParser
from itertools import izip
from qiime.golay import decode_golay_12, get_invalid_golay_barcodes
from qiime.parse import parse_mapping_file_to_dict
from qiime.split_libraries import check_map, expand_degeneracies
from qiime.split_libraries_fastq import correct_barcode, FastqParseError
from skbio.parse.sequences import parse_fastq
from qiime.util import create_dir, qiime_system_call, get_qiime_temp_dir
import re
import tempfile
import os


class PairedEndParseError(Exception):
    pass


class PrimerMismatchError(Exception):
    pass


class LowConsensusScoreError(Exception):
    pass


class SeqLengthMismatchError(Exception):
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
        raise PrimerMismatchError

    before_primer = seq[:primer_idx]
    after_primer = seq.replace(before_primer + primer, '', 1)
    return before_primer, primer, after_primer


def read_input_file(sequence_read_fps, mapping_fp, output_dir,
                    barcode_type, max_barcode_errors,
                    min_consensus, max_cluster_ratio, option_parser):
    """
    Reads mapping file, input file, and other command line arguments
    fills dictionary called consensus_seq_lookup which will contain:
    sample ID -> random barcode -> consensus_seq
    """

    random_bc_lookup = defaultdict(lambda:
                                   defaultdict(lambda:
                                               defaultdict(int)))

    consensus_seq_lookup = defaultdict(lambda:
                                       defaultdict(lambda:
                                                   defaultdict(int)))

    BARCODE_COLUMN = 'BarcodeSequence'
    REVERSE_PRIMER_COLUMN = 'ReversePrimer'

    if barcode_type == 'golay_12':
        barcode_correction_fn = decode_golay_12
        barcode_len = 12
    else:
        barcode_correction_fn = None

        try:
            barcode_len = int(barcode_type)
        except ValueError:
            option_parser.error("Invalid barcode type '%s'. The barcode type "
                                "must be either golay_12 or a positive "
                                "integer indicating the barcode length." %
                                barcode_type)

    if max_barcode_errors < 0:
        option_parser.error("--max_barcode_errors must be greater than or "
                            "equal to zero. You provided %.4f." %
                            max_barcode_errors)

    if barcode_len < 1:
        option_parser.error("Invalid barcode length: %d. Must be greater "
                            "than zero." % barcode_len)

    seq_fps = sequence_read_fps

    if len(seq_fps) != 2:
        option_parser.error("You must provide exactly two sequence read "
                            "filepaths, the first for forward reads and "
                            "second for reverse reads. You specified %d "
                            "filepaths." % len(seq_fps))

    with open(mapping_fp, 'U') as map_f:
        #  Ensures that sample IDs and barcodes are unique, that barcodes are
        #  all the same length, and that primers are present. Ensures barcodes
        #  and primers only contain valid characters.
        _, _, bc_to_sid, _, _, bc_to_fwd_primers, _ = check_map(map_f, False)
        map_f.seek(0)

        #  TODO: add reverse primer validation similar to what check_map does
        # (probably just modify check_map to account for reverse primer).

        metadata_map = parse_mapping_file_to_dict(map_f)[0]
        bc_to_rev_primers = {}
        for sid, md in metadata_map.items():
            if REVERSE_PRIMER_COLUMN in md:
                bc_to_rev_primers[md[BARCODE_COLUMN]] = expand_degeneracies
                (md[REVERSE_PRIMER_COLUMN].upper().split(','))
            else:
                option_parser.error("The %s column does not exist in the "
                                    "mapping file. %s is required." %
                                    (REVERSE_PRIMER_COLUMN,
                                     REVERSE_PRIMER_COLUMN))

    #  Make sure our barcodes(which are guaranteed to be the same length at
    #  this point) are the correct length that the user specified.
    barcode_len_in_map = len(bc_to_sid.keys()[0])
    if barcode_len_in_map != barcode_len:
        option_parser.error("Barcodes in mapping file are of length %d, but "
                            "expected barcodes of length %d." %
                            (barcode_len_in_map, barcode_len))

    if barcode_type == 'golay_12':
        invalid_golay_barcodes = get_invalid_golay_barcodes(bc_to_sid.keys())

        if invalid_golay_barcodes:
            option_parser.error("Some or all barcodes in the mapping file are "
                                "not valid golay codes. Do they need to be "
                                "reverse complemented? If these are not golay "
                                "barcodes pass --barcode_type 12 to disable "
                                "barcode error correction, or pass "
                                "--barcode_type  # if the barcodes are not 12 "
                                "base pairs, where   #  is the size of the "
                                "barcodes.\n\nInvalid barcodes: %s" %
                                ' '.join(invalid_golay_barcodes))

    header_idx = 0
    seq_idx = 1
    qual_idx = 2
    barcode_errors_exceed_max_count = 0
    barcode_not_in_map_count = 0
    primer_mismatch_count = 0
    fwd_read_f = open(seq_fps[0], 'U')
    rev_read_f = open(seq_fps[1], 'U')
    random_bcs = list()
    for fwd_read, rev_read in izip(
            parse_fastq(fwd_read_f, strict=False),
            parse_fastq(rev_read_f, strict=False)):
            # Confirm match between read headers.

        if fwd_read[header_idx] != rev_read[header_idx]:
            raise PairedEndParseError("Headers of forward and reverse reads "
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
            continue  # not sure about indentation

        # Extract the random barcode and primer from the forward read.
        possible_primers = bc_to_fwd_primers[corrected_barcode].keys()

        try:
            random_bc, _, clean_fwd_seq = extract_primer(fwd_seq,
                                                         possible_primers,
                                                         min_idx=5,
                                                         max_idx=20)
            random_bcs.append(random_bc)
        except PrimerMismatchError:
            primer_mismatch_count += 1
            continue

        possible_primers = bc_to_rev_primers[corrected_barcode]

        try:
            # phase_seq, _, clean_rev_seq = extract_primer(rev_seq,
            # possible_primers)
            pass
        except PrimerMismatchError:
            primer_mismatch_count += 1
            # continue
        clean_rev_seq = fwd_seq

        random_bc_lookup[sample_id][random_bc][(clean_fwd_seq, clean_rev_seq)] += 1
    cluster_ratios = list()
    temp_dir = get_qiime_temp_dir
    random_bc_keep = select_unique_rand_bcs(random_bcs)
    for sample_id in random_bc_lookup:
        for random_bc in random_bc_lookup[sample_id]:
            if random_bc_keep[random_bc] == 1:
                fasta_tempfile = tempfile.NamedTemporaryFile(dir=output_dir, delete=False)
                fasta_tempfile_name = fasta_tempfile.name
                max_freq = 0
                for seq_count_this_barcode, fwd_rev_seq in enumerate(random_bc_lookup[sample_id][random_bc]):
                    fwd_seq, rev_seq = fwd_rev_seq
                    p_line = ">" + random_bc + "|" + str(random_bc_lookup[sample_id][random_bc][fwd_rev_seq]) + "\n" + fwd_seq + "\n"
                    fasta_tempfile.write(p_line)
                    num_seq_this_barcode = seq_count_this_barcode
                    if random_bc_lookup[sample_id][random_bc][fwd_rev_seq] > max_freq:
                        max_freq = random_bc_lookup[sample_id][random_bc][fwd_rev_seq]
                        majority_seq = fwd_seq
                fasta_tempfile = open(fasta_tempfile_name, 'r')
                cluster_ratio = get_cluster_ratio(fasta_tempfile_name)
                # name is passed because there is system call inside the function
                # function does not open the file
                cluster_ratios.append(cluster_ratio)
                if cluster_ratio < max_cluster_ratio:
                    consensus_seq = majority_seq
                else:
                    consensus_seq = get_consensus(fasta_tempfile, min_consensus)
                consensus_seq_lookup[sample_id][random_bc] = consensus_seq
                fasta_tempfile.close()
                os.unlink(fasta_tempfile_name)

    fwd_read_f.close()
    rev_read_f.close()
    return consensus_seq_lookup


def get_cluster_ratio(temp_file):
    """
    Uses uclust to calculate cluster ratio
    cluster_ratio=num_of_seq_in_cluster_with_max_seq/num_of_seq_in cluster_with_second_higest_seq
    """
    temp_dir = get_qiime_temp_dir()
    uclust_tempfile = tempfile.NamedTemporaryFile(dir=temp_dir, mode='w', delete=False)
    uclust_tempfile_name = uclust_tempfile.name
    command = "uclust --usersort --input " + temp_file +\
              " --uc " + uclust_tempfile_name + " --id 0.98 --log log"
    qiime_system_call(command)
    uclust_tempfile.close()
    uclust_tempfile = open (uclust_tempfile_name, 'r')
    count = 0
    count_lookup = {}

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

    sorted_counts_in_clusters = sorted(count_lookup.iteritems(), key=lambda x: x[1])
    uclust_tempfile.close()    
    os.unlink(uclust_tempfile.name)
    try:
        return float(str(sorted_counts_in_clusters[0][1]))/float(str(sorted_counts_in_clusters[1][1]))
    except IndexError:
        return 1


def get_consensus(fasta_tempfile, min_consensus):
    """
    Returns consensus sequence from a set of sequences
    input: fasta file, min_consensus
    fasta_file should be in the following format:
    >random_bc|number
    seq
    >random_bc|number
    seq
    ....

    where: number = number of times the particular seq has appeared with this random_barcode
    """

    length = 0
    seqs = []  # list of sequences in fasta file
    counts = []  # count of occurence of sequnces
    this_seq_count = -1  # counter for sequence
    number_of_seqs = 0
    # I have not used a function from skitbio for reading the
    # fasta file, because of the particular format of fasta file.
    # i.e. random barcode, and count in ID line
    
    fasta_tempfile_name =fasta_tempfile.name

    # read the fasta file
    # store seqs and counts
    for line in fasta_tempfile:
        if re.search(r'\>', line):
            this_seq_count += 1
            seqs.append("")
            counts.append(0)
            RE_output = re.search(r'\>\w+\|(\d+)', line)
            counts[this_seq_count] = RE_output.group(1)
            counts[this_seq_count] = int(counts[this_seq_count])
            number_of_seqs += counts[this_seq_count]
        else:
            line = line.rstrip('\n')
            try:
                seqs[this_seq_count] = seqs[this_seq_count] + line
            except IndexError:
                seqs[this_seq_count] = line

    length = len(seqs[0])
    # if length != old_length:
    # print length, old_length
    # raise SeqLengthMismatchError
    # lookup is a 2D list
    # It has N rows and M columns
    # N: number of seqs
    # M: length of the seqs(should be same)

    lookup = {}
    for i in range(length):
        lookup[i] = {}
        for j in range(this_seq_count):
            lookup[i][j] = 0

    for i in range(length):
        for this_seq_count, seq in enumerate(seqs):
            try:
                lookup[i][seq[i]] += counts[this_seq_count]
            except KeyError:
                lookup[i][seq[i]] = counts[this_seq_count]

    consensus = ''      # consesus sequence
    con_score = ''
    # consensus score: string. At each position: range 0-9
    # at each position, 10 * occurence of max base / number_of_seqs
    # 10 is converted to 9

    count = 0

    for index in range(length):
        sorted_bases = sorted(lookup[index].iteritems(), key=lambda x: x[1])
        max_base = str(sorted_bases[-1])

        (base, num) = max_base.split(',')
        num = num.replace(")", "")
        base = base.replace("(", "")
        base = base.replace("'", "")
        num = int(num)
        score = float(10 * float(num) / number_of_seqs)

        if score >= 9:
            score = 9

        if score < min_consensus:
            raise LowConsensusScoreError

        score = str(int(score))
        consensus = consensus+base
        con_score += score
        count += 1

    if con_score >= min_consensus:
            return consensus


def select_unique_rand_bcs(rand_bcs):
    """
    Attempts to select true barcodes from set of barcodes
    i.e. removes barcodes that might be artifacts
    due to sequencing errors.
    Uses uclust to remove barcodes that are similar thatn
    threshold.
    returns: a dictionary which has value of 1 for unique
    barcodes, and 0 for others.
    """
    unique_threshold = 0.86
    unique_rand_bcs = {}
    temp_dir = get_qiime_temp_dir()

    fasta_tempfile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, mode='w')
    uclust_tempfile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, mode='r')
    uclust_tempfile_name = uclust_tempfile.name
    fasta_tempfile_name = fasta_tempfile.name

    p_line = ""
    for count_rand_bc, rand_bc in enumerate(rand_bcs):
        p_line = p_line + ">" + str(rand_bc) + "\n" + rand_bc + "\n"
        unique_rand_bcs[rand_bc] = 0
    fasta_tempfile.write(p_line)
    qiime_system_call("uclust --usersort --input " + fasta_tempfile.name
                      + " --uc " + uclust_tempfile.name + " --id " +
                      str(unique_threshold) + " --log log")
    for line in uclust_tempfile:
        if re.search('^C', line):
            pieces = line.split('\t')
            unique_rand_bc = pieces[8]
            unique_rand_bcs[unique_rand_bc] = 1
    fasta_tempfile.close()
    uclust_tempfile.close()
    os.unlink(fasta_tempfile_name)
    os.unlink(uclust_tempfile_name)

    return unique_rand_bcs
