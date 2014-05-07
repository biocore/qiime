from collections import defaultdict
from itertools import izip
from qiime.golay import decode_golay_12, get_invalid_golay_barcodes
from qiime.parse import  parse_mapping_file_to_dict
from qiime.split_libraries import check_map, expand_degeneracies
from qiime.split_libraries_fastq import correct_barcode, FastqParseError
from qiime.util import create_dir, qiime_system_call, get_qiime_temp_dir
import re
import tempfile
import skbio.parse.sequences.parse_fastq

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

            if (min_idx is not None and primer_idx < min_idx) or \
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
                    barcode_type, max_barcode_errors, min_consensus, max_cluster_ratio):
    """
    Reads mapping file, input file, and other command line arguments
    fills dictionary called consensus_seq_lookup which will contain:
    sample ID -> random barcode -> consensus_seq
    """

    random_bc_lookup = defaultdict(lambda:
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

    create_dir(output_dir)

    with open(mapping_fp, 'U') as map_f:
        #  Ensures that sample IDs and barcodes are unique, that barcodes are
        #  all the same length, and that primers are present. Ensures barcodes
        #  and primers only contain valid characters.
        _, _, bc_to_sid, _, _, bc_to_fwd_primers, _ = check_map(map_f, False)
        map_f.seek(0)

        #  TODO: add reverse primer validation similar to what check_map does
        #  (probably just modify check_map to account for reverse primer).
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

    #  Make sure our barcodes (which are guaranteed to be the same length at
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

    for fwd_read, rev_read in izip(
            skbio.parse.sequences.parse_fastq(fwd_read_f, strict=False),
            skbio.parse.sequences.parse_fastq(rev_read_f, strict=False)):
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

        #  Correct the barcode (if applicable) and map to sample ID.
        num_barcode_errors, corrected_barcode, _, sample_id = correct_barcode\
		(barcode, bc_to_sid, barcode_correction_fn)

        #  Skip barcodes with too many errors.
        if num_barcode_errors > max_barcode_errors:
            barcode_errors_exceed_max_count += 1
            continue


        # Skip unassignable reads.
        # TODO: do we want to keep unassignable reads?
        # If so, how to choose the primer to remove?


        if sample_id is None:
            barcode_not_in_map_count += 1
            continue  # not sure about indentation

        # TODO: other quality filtering using
        # qiime.split_libraries_fastq.quality_filter_sequence
        # Extract the random barcode and primer from the forward read.
        possible_primers = bc_to_fwd_primers[corrected_barcode].keys()

        try:
            random_bc, _, clean_fwd_seq = extract_primer(fwd_seq,
                                                         possible_primers)

        except Exception:   # PrimerMismatchError:
            primer_mismatch_count += 1
            continue

        possible_primers = bc_to_rev_primers[corrected_barcode]

        try:
            phase_seq, _, clean_rev_seq = extract_primer(rev_seq,
                                                         possible_primers)

        except Exception:  # PrimerMismatchError
			primer_mismatch_count += 1
			continue
		



		random_bc_lookup[sample_id][random_bc][(clean_fwd_seq, clean_rev_seq)] += 1

	temp_dir=get_qiime_temp_dir
	fasta_tempfile = tempfile.NamedTemporaryFile()
	for sample_id in random_bc_lookup:
		for random_bc in random_bc_lookup[sample_id]:
			max_freq=0
			for seq_count_this_barcode, fwd_rev_seq in enumerate(random_bc_lookup[sample_id][random_bc]):
				fwd_seq, rev_seq = fwd_rev_seq
                fasta_tempfile.write(">" + random_bc_lookup[sample_id][random_bc][fwd_rev_seq] + "\n" + fwd_seq+"\n")
				num_seq_this_barcode=seq_count_this_barcode
				if random_bc_lookup[sample_id][random_bc][fwd_rev_seq] > max_freq:
					max_freq=random_bc_lookup[sample_id][random_bc][fwd_rev_seq]
					majority_seq=fwd_seq
			cluster_ratio = get_cluster_ratio(fasta_tempfile.name)  
			# name is passed instead of file handle because
			# there is a system call inside the function

			if cluster_ratio < max_cluster_ratio:
				consensus_seq = majority_seq
			else:
				consensus_seq = get_consensus(fasta_tempfile, min_consensus) # read append mode?  
			consensus_seq_lookup[sample_id][random_bc]=consensus_seq

    fwd_read_f.close()
    rev_read_f.close()
    return consensus_seq_lookup

    # if the ratio is higher than threshold,
    # we need to create a consensus sequence for the given sequence.
    # threshold= 2.5


def get_cluster_ratio(temp_file):
    """
    uses uclust to calculate cluster ratio	
	cluster_ratio=num_of_seq_in_cluster_with_max_seq/num_of_seq_in cluster_with_second_higest_seq
    """

    uclust_tempfile = tempfile.NamedTemporaryFile()
    qiime_system_call
    ("uclust --usersort --input " + temp_file.name
     + " --uc "+uclust_tempfile.name + "--id 0.3 --log log")
    count = 0
    count_lookup = {}

    for line in uclust_tempfile.name:
        if re.search("^S|^H", line):
            pieces = line.split('\t')
            try:
                count_lookup[pieces[1]] += 1
            except KeyError:
                count_lookup[pieces[1]] = 1
            count += 1

    sorted_counts_in_cluters = sorted(count_lookup.iteritems(), key=lambda x: x[1])
    try:
        sorted_counts_in_cluters[-2]
    except IndexError:
        return 1

    return sorted_counts_in_cluters[-1]/sorted_counts_in_cluters[-2]


def get_consensus(fasta_tempfile, min_consensus):
    """
    returns consensus sequence from a set of sequences
	input: fasta file, min_consensus
	fasta_file should contain count of the particular sequence
	after the > sign
	
    """
	length = 0
	seqs = list()
	counts = list()
	seq = ""
	# read the fasta file
	# append all sequences in a list called seqs
	# count of each seq is saved in counts

	for line in fasta_tempfile:
		if re.serch('\>', line):
			if not seq is "":
				old_length = length
				length = len(seq)
				if length != old_length:
					raise SeqLengthMismatchError 
							
				seqs.append(seq)			
				counts.append(count)
			RE_output=re.serch('\>(\d+)', line)
			count=RE_output.group('1')
			seq=""
		else:    
			seq = seq+line
			seq = seq.rstrip('\n')
	
	lookup = {}
	for i in range(length):
        lookup[i] = {}

    for i in range(length):
        for seq in seqs:
            try:
                lookup[i][seq[i]]
            except KeyError:
                lookup[i][seq[i]] = 1
            else:
                lookup[i][seq[i]] += 1

    consensus = ''      # consesus sequence
    con_score = ''
    # consensus score: string. At each position: range 1-10.
    # at each position, 10 * occurence of max base / number_of_seqs
    # 10 is converted to 0
    count = 0

    for index in lookup:
        sorted_bases = sorted(lookup[index].iteritems(), key=lambda x: x[1])
        max_base = str(sorted_bases[-1])

        (base, num) = max_base.split(',')
        num = num.replace(")", "")
        base = base.replace("(", "")
        base = base.replace("'", "")
        num = int(num)
        score = 10 * num / number_of_seqs

        if score < 1:
            score = 1

        if score == 10:
            score = 0

        score = str(score)
        consensus = consensus+base
        con_score += score
        count += 1
	if con_score >= min_consensus:
    		return consensus
	else:
		raise LowConsensusScoreError


def select_unique_rand_bcs(rand_bcs):
	"""
	attempts to select true barcodes from set of barcodes
	i.e. removes barcodes that might be artifacts
	due to sequencing errors.
	"""

	unique_rand_bcs=list()
	temp_dir=get_qiime_temp_dir
	fasta_tempfile = tempfile.NamedTemporaryFile(dir=temp_dir)
	uclust_tempfile = tempfile.NamedTemporaryFile(dir=temp_dir)	
	for count_rand_bc, rand_bc in enumerate(rand_bcs):
		fasta_tempfile.write(">" + str(count_rand_bc) +"\n"+rand_bc+"\n")
	
	qiime_system_call("uclust --usersort --input " + fasta_tempfile.name + " --uc "+uclust_tempfile.name + "--id 0.3 --log log")	
	for line in uclust_tempfile.name:
		if re.search('^C', line):
			pieces = line.split('\t')
			unique_rand_bc, num_of_reads = pieces[8].split('\.') 
			unique_rand_bcs.append(uniue_rand_bc)
	fasta_tempfile.close()
	uclust_tempfile.close()
	return unique_rand_bcs

class PairedEndParseError:
    pass

class PrimerMismatchError(Exception):
    pass

class LowConsensusScoreError(Exception):
    pass

class SeqLengthMismatchError(Exception):
	pass
