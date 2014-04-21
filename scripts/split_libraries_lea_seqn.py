#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout", "Charudatta Navare"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Charudatta Navare"
__email__ = "charudatta.navare@gmail.com"

from collections import defaultdict
from itertools import izip
from qcli import parse_command_line_parameters, make_option
from qiime.golay import decode_golay_12, get_invalid_golay_barcodes
from qiime.parse import MinimalFastqParser, parse_mapping_file_to_dict
from qiime.split_libraries import check_map, expand_degeneracies
from qiime.split_libraries_fastq import correct_barcode, FastqParseError
from qiime.util import create_dir, qiime_system_call
from lib_split_libraries_lea_seq import get_cluster_ratio, consensus_seq, extract_primer, select_majority_sequence, get_consensus, read_input_file
import tempfile
import re

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
	make_option('--consensus_threshold', type='float',
                help='threshold for consensus score'
		'the minimum score allowable at any position in sequence'
		'where the score is calulated as:'
		'occurence of base in consensus sequence/ total sequences'
                     '[default: %default]',
                default=0.66)
	make_option('--threshold_for_cluster_ratio', type='float',
                help='threshold for cluster ratio'
			'the maximum allowable cluster ratio'
			'above which you need to find the consensus sequence'
			'for the given sequences'
                     '[default: %default]',
                default=2.5)

]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    barcode_type = opts.barcode_type
    max_barcode_errors = opts.max_barcode_errors
    output_dir=opts.output_dir
    mapping_fp=opts.mapping_fp
    sequence_read_fps=opts.sequence_read_fps
    consensus_threshold=opts.consensus_threshold
    cluster_ratio_threshold=opts.threshold_for_cluster_ratio

    random_bc_lookup = defaultdict(lambda:
                                   defaultdict(lambda:
                                               defaultdict(int)))

    read_input_file(sequence_read_fps, mapping_fp, output_dir, barcode_type, max_barcode_errors, random_bc_lookup)
    consensus_outfile = open(os.path.join(output_dir, mapping_fp+"_consensus.fasta"), "w")

    bc_count=1	# counter for random barcode
	
    for sample_id in random_bc_lookup:
	for random_bc in random_bc_lookup[sample_id]:

		bc_count += 1 
		seq_count_this_barcode = 0	#counter for sequences associated with this random barcode
		fasta_tempfile = tempfile.NamedTemporaryFile()
		for fwd_rev_seq in random_bc_lookup[sample_id][random_bc]:

			seq_count_this_barcode += 1

			str_seq = str(fwd_rev_seq)	
			#str_seq has forward and reverse sequence in string format, fwd and rev seqparated by comma
			(fwd_seq, rev_seq) = str_seq.split(",")
			
			fwd_seq = fwd_seq.replace("(","")
	    		fwd_seq = fwd_seq.replace("'","")
			rev_seq = rev_seq.replace("(","")
	    		rev_seq = rev_seq.replace("'","")
			
			fasta_tempfile.write(">"+str(seq_count_this_barcode)+"\n"+fwd_seq+"\n") # write sequences in fasta format, 

		cluster_ratio = get_cluster_ratio(fasta_tempfile.name)
		fasta_seqs = get_seqs(fasta_tempfile.name)

		if cluster_ratio < cluster_ratio_threshold:
			consensus_seq=select_majority_sequence(fasta_seqs)
		else:
			consensus_seq=get_consensus_seq(fasta_seqs, consensus_threshold)

		consensus_outfile.write(">"+str(bc_count)+"\n"+consensus_seq+"\n")	

if __name__ == "__main__":
    main()
