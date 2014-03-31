#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout", "Charudatta Navare"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from collections import defaultdict
from itertools import izip
from qcli import parse_command_line_parameters, make_option
from qiime.golay import decode_golay_12, get_invalid_golay_barcodes
from qiime.parse import MinimalFastqParser, parse_mapping_file_to_dict
from qiime.split_libraries import check_map, expand_degeneracies
from qiime.split_libraries_fastq import correct_barcode, FastqParseError
from qiime.util import create_dir
import os
import re


#########################
def cluster_ratio(seqs):
"""
cluster ratio using uclust: if the ratio is high,
we need to create a consensus sequence for the given sequence.
threshold= 2.5

"""

	seq_str = "";
	seq_count=0;

	for s in seqs:
		seq_str = seq_str+">"+ str(seq_count) + "."+ str(seq_count)+"\n"+s+"\n";
		seq_count+=1

	seq_out = "sequences.fas";
	seq_file = open (seq_out, "w");
	seq_file.write(seq_str)

	uc_out = "seq.uc";
	os.system("uclust --usersort --input ucGLOBAL_COUNT.fa --uc seq.uc --id 0.3 --log log");
	
	count=0;
	count_hash={}
	clust_infile=open (uc_out, "r");
	for line in clust_infile:
		if re.search("^S|^H", line):
			
			pieces = line.split('\t');
			(num_reads, id) = pieces[8].split('.');
			count_hash[pieces[1]]+=num_reads;
			count+=1

	key2=sorted(count_hash.iteritems(), key=lambda x:x[1])
	
	return key2


def consensus_seq(seqs):
"""
returns consensus sequence for the given sequences

"""

    total=len(seqs)

	#string=len(seqs[0])
	#print type(string)


	length=len(seqs[0])
	lookup={}
	for i in range(length):
	    lookup[i]={}

	#lookup={0:{}, 1: {}, 2:{}, 3:{}, 4:{}, 5:{}}

	for i in range(length):
		for seq in seqs:
			try:
				lookup[i][seq[i]]
			except KeyError:
				lookup[i][seq[i]]=1
			else:	
				lookup[i][seq[i]]+=1
	consensus=''
	con_score=''
	count=0
	for key1 in lookup:
	    key2=sorted(lookup[key1].iteritems(), key=lambda x:x[1])
	    #print key1, key2[-1]
	    key3=str(key2[-1])

	    (base, num)=key3.split(',')
	    num=num.replace(")","")
	    base=base.replace("(","")
	    base=base.replace("'","")

	    num=int(num)
	    score= 10* num /total
	    if score ==10:
		score =0
	    #    print score
	    score=str(score)
	    consensus=consensus+base
	    con_score+=score
	    count+=1


	return consensus, con_score

def extract_primer(seq, possible_primers, min_idx=None, max_idx=None):
    """

``min_idx`` and ``max_idx`` are inclusive.

TODO: allow primer mismatches?
TODO: this may be dangerous as multiple primers may match, and we're
passed these in random order.
"""
    primer_idx = None
    primer = None
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
        # TODO: fill in error message
        raise PrimerMismatchError

    before_primer = seq[:primer_idx]

    # TODO: what to do with random barcodes/phases that have ambiguous bases?
    # TODO: seq.replace may be dangerous here- should verify it won't remove
    # the cruft from somewhere unintended.
    return before_primer, primer, seq.replace(before_primer + primer, '', 1)


###########
script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(("","",""))
script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-i', '--sequence_read_fps', type='existing_filepaths',
                help='the forward and reverse sequence read fastq files '
                     '(comma-separated)'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='directory to store output files'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='metadata mapping file')
]
script_info['optional_options'] = [
    make_option('--barcode_type', type='string',
                help='the type of barcode used. This can be an integer, e.g. '
                     '6 for length 6 barcodes, or golay_12 for golay error-'
                     'correcting barcodes. Error correction will only be '
                     'applied for golay_12 barcodes [default: %default]',
                default='golay_12'),
    make_option('--max_barcode_errors', type='float',
                help='the maximum allowable number of errors in the barcode '
                     'if passing --barcode_type golay_12 [default: %default]',
                default=1.5)
]
script_info['version'] = __version__

BARCODE_COLUMN = 'BarcodeSequence'
REVERSE_PRIMER_COLUMN = 'ReversePrimer'

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    barcode_type = opts.barcode_type
    max_barcode_errors = opts.max_barcode_errors

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

    seq_fps = opts.sequence_read_fps

    if len(seq_fps) != 2:
        option_parser.error("You must provide exactly two sequence read "
                            "filepaths, the first for forward reads and "
                            "second for reverse reads. You specified %d "
                            "filepaths." % len(seq_fps))

    create_dir(opts.output_dir)

    with open(opts.mapping_fp, 'U') as map_f:
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
    fwd_read_f = open(seq_fps[0], 'U')
    rev_read_f = open(seq_fps[1], 'U')

    barcode_errors_exceed_max_count = 0
    barcode_not_in_map_count = 0
    primer_mismatch_count = 0

    # sample ID -> random barcode -> (fwd seq, rev seq) -> count
    random_bc_lookup = defaultdict(lambda:
                                   defaultdict(lambda:
                                               defaultdict(int)))

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

        # Skip unassignable reads. TODO: do we want to keep unassignable reads?
        # If so, how to choose the primer to remove?
        if sample_id is None:
          barcode_not_in_map_count += 1

          continue

        # TODO: other quality filtering using
        # qiime.split_libraries_fastq.quality_filter_sequence

        # Extract the random barcode and primer from the forward read.
        possible_primers = bc_to_fwd_primers[corrected_barcode].keys()

        # TODO: allow user to parameterize the min and max random barcode
        # lengths.
        try:
            random_bc, _, clean_fwd_seq = extract_primer(fwd_seq,
                                                         possible_primers,
                                                         min_idx=16,
                                                         max_idx=18)

	except PrimerMismatchError:
            primer_mismatch_count += 1
            continue

        # TODO: how to handle phasing? Should the user be responsible for
        # providing a standard sequence length that we truncate to here?

        # Clean up reverse read by extracting the phased portion and the
        # reverse primer.
        possible_primers = bc_to_rev_primers[corrected_barcode]

        # TODO: allow user to parameterize the phasing length.
        try:
            phase_seq, _, clean_rev_seq = extract_primer(rev_seq,
                                                         possible_primers,
                                                         max_idx=3)

        except PrimerMismatchError:
            # TODO: count fwd and rev mismatches differently?
            primer_mismatch_count += 1

            continue

        random_bc_lookup[sample_id][random_bc][(clean_fwd_seq, clean_rev_seq)] += 1

    fwd_read_f.close()
    rev_read_f.close()

    for (k,v) in random_bc_lookup:
        print k, ":\n"


    print barcode_errors_exceed_max_count
    print barcode_not_in_map_count
    print primer_mismatch_count
    print
    print random_bc_lookup



# use uclust to get clusters
# if the clusters are more, indicating low quality, use consensus sequence and majority sequence
# if the clusrers are low, (high quality), select the best sequence as consensus majority


# print seqs in txt file
# run uclust at 97%
# cluster ratio
# run uclust to get unique barcodes
# only select unique barcodes


class PairedEndParseError(FastqParseError):
    pass

class PrimerMismatchError(Exception):
    pass



if __name__ == "__main__":
    main()

