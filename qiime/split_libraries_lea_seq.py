from collections import defaultdict
from itertools import izip
from qcli import parse_command_line_parameters, make_option
from qiime.golay import decode_golay_12, get_invalid_golay_barcodes
from qiime.parse import MinimalFastqParser, parse_mapping_file_to_dict
from qiime.split_libraries import check_map, expand_degeneracies
from qiime.split_libraries_fastq import correct_barcode, FastqParseError
from qiime.util import create_dir, qiime_system_call
import re
import tempfile

def extract_primer(seq, possible_primers, min_idx=None, max_idx=None):
    """
	Extracts primers from sequence, given possible primers
	returns before_primer, primer, sequence without primer

    """

    primer_idx = None #index of primer in sequence
    primer = None	#primer in sequence

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
    pruned_sequence=seq.replace(before_primer + primer, '', 1)
    return before_primer, primer, pruned_sequence 



def read_input_file(sequence_read_fps, mapping_fp, output_dir, barcode_type, max_barcode_errors, random_bc_lookup):
    '''
	Reads mapping file, input file, and other command line arguments
	fills dictionary called random_bc_lookup which will contain: 
	sample ID -> random barcode -> (fwd seq, rev seq) -> count
    '''

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
        # Ensures that sample IDs and barcodes are unique, that barcodes are
        # all the same length, and that primers are present. Ensures barcodes
        # and primers only contain valid characters.
        _, _, bc_to_sid, _, _, bc_to_fwd_primers, _ = check_map(map_f, False)
        map_f.seek(0)

        # TODO: add reverse primer validation similar to what check_map does
        # (probably just modify check_map to account for reverse primer).
        metadata_map = parse_mapping_file_to_dict(map_f)[0]
        bc_to_rev_primers = {}
        for sid, md in metadata_map.items():
            if REVERSE_PRIMER_COLUMN in md:
                bc_to_rev_primers[md[BARCODE_COLUMN]] = expand_degeneracies(
                        md[REVERSE_PRIMER_COLUMN].upper().split(','))
            else:
                option_parser.error("The %s column does not exist in the "
                                    "mapping file. %s is required." %
                                    (REVERSE_PRIMER_COLUMN,
                                     REVERSE_PRIMER_COLUMN))

    # Make sure our barcodes (which are guaranteed to be the same length at
    # this point) are the correct length that the user specified.
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
                                "--barcode_type # if the barcodes are not 12 "
                                "base pairs, where # is the size of the "
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
            MinimalFastqParser(fwd_read_f, strict=False),
            MinimalFastqParser(rev_read_f, strict=False)):
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

        # Grab the barcode sequence. It is always at the very end of the
        # forward read. Strip the barcode from the sequence.
        barcode = fwd_seq[-barcode_len:]
        fwd_seq = fwd_seq[:-barcode_len]

        # Correct the barcode (if applicable) and map to sample ID.
        num_barcode_errors, corrected_barcode, _, sample_id = correct_barcode(
                barcode, bc_to_sid, barcode_correction_fn)

        # Skip barcodes with too many errors.
        if num_barcode_errors > max_barcode_errors:
          barcode_errors_exceed_max_count += 1
          continue

        # Skip unassignable reads. 
        if sample_id is None:
          barcode_not_in_map_count += 1

          continue

        # Extract the random barcode and primer from the forward read.
        possible_primers = bc_to_fwd_primers[corrected_barcode].keys()

        try:
            random_bc, _, clean_fwd_seq = extract_primer(fwd_seq,
                                                         possible_primers,
                                                         min_idx=16,
                                                         max_idx=18)

	except Exception: #PrimerMismatchError:
            primer_mismatch_count += 1
            continue

        possible_primers = bc_to_rev_primers[corrected_barcode]

        try:
            phase_seq, _, clean_rev_seq = extract_primer(rev_seq,
                                                         possible_primers)

        except Exception:#PrimerMismatchError
            primer_mismatch_count += 1

            continue

        random_bc_lookup[sample_id][random_bc][(clean_fwd_seq, clean_rev_seq)] += 1
        #sample ID -> random barcode -> (fwd seq, rev seq) -> count

	fwd_read_f.close()
	rev_read_f.close()
	return 0	    



def get_cluster_ratio(temp_file):
	"""
		cluster ratio=
		number of sequences in minority cluster/number of sequences in majority cluster  
		uclust is used to cluster the sequences 
		input: fasta file
	"""

	uclust_tempfile=tempfile.NamedTemporaryFile()
	qiime_system_call("uclust --usersort --input "+temp_file.name+" --uc "+uclust_tempfile.name+"--id 0.3 --log log");
	
	count=0;
	seqs_in_cluster={}

	for line in uclust_tempfile.name:
		if re.search("^S|^H", line):
			
			pieces = line.split('\t');
			try:			
				seqs_in_cluster[pieces[1]]+=1;
			except KeyError:
				seqs_in_cluster[pieces[1]]=1;
			count+=1

	sorted_seqs_in_cluster = sorted(seqs_in_cluster.iteritems(), key=lambda x:x[1])
	try :
		sorted_seqs_in_cluster[1]		
	except IndexError:
		return 1

	return sorted_seqs_in_cluster[0]/sorted_seqs_in_cluster[1]


def get_consensus_seq(seqs, consensus_threshold):
	"""
	returns consensus sequence from a set of sequences
	input: a list of sequences, and the minimum allowable consensus score	
	"""

	number_of_seqs=len(seqs)
	length=len(seqs[0])	#assumes that sequences will have same length 

	lookup={}
	for i in range(length):
	    lookup[i]={}

	for i in range(length):
		for seq in seqs:
			try:
				lookup[i][seq[i]]
			except KeyError:
				lookup[i][seq[i]]=1
			else:	
				lookup[i][seq[i]]+=1


	consensus=''	#consesus sequence
	con_score=''	#consensus score: string depicting consensus score at each position. range 1-10.
			# at each position, 10 * occurence of base that occured max times / number_of_seqs
			# 10 is converted to 0
	count=0

	for index in lookup:
	    sorted_bases=sorted(lookup[index].iteritems(), key=lambda x:x[1])
	    max_base=str(sorted_bases[-1])

	    (max_base, max_num)=max_base.split(',')
	    max_num=max_num.replace(")","")
	    max_base=max_base.replace("(","")
	    max_base=base.replace("'","")

	    num=math.ceil(num)
	    score= 10* num /number_of_seqs

	    if score ==10:
		score =0

	    score=str(score)
	    consensus=consensus+base
	    con_score+=score
	    count+=1
	return consensus, con_score


def get_seqs(fasta_file_name):
	#return seqs from fasta file	
	pass

def select_majority_sequence(seqs):
	#return majority seq from set of sequence
	pass

def select_unique_primers:
	#return unique primers(which are different from each other, threshold=)
	pass

class PairedEndParseError:
    pass

class PrimerMismatchError(Exception):
    pass
