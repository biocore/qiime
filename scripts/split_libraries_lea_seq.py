#!/usr/bin/env python
from __future__ import division

__author__ = "Charudatta Navare"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Charudatta Navare", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Charudatta Navare"
__email__ = "charudatta.navare@gmail.com"

from qcli import parse_command_line_parameters, make_option
from qiime.split_libraries_lea_seq import (get_cluster_ratio, select_majority_sequence, get_consensus, read_input_file)
import tempfile

script_info={}
script_info['brief_description'] = "Implements Low-Error Amplicon Sequencing (LEA-Seq)"

script_info['script_description'] = """Implements Low-Error Amplicon Sequencing (LEA-Seq) method, \
described in: Faith, Jeremiah J., et al. \
The long-term stability of the human gut microbiota. Science 341.6141 (2013).
This method is based on redundant sequencing of a set of linear PCR template extensions of 16S rRNA genes. The oligonucleotide primer that is used for PCR template extensions is labeled with a random barcode 5' to the universal 16S rRNA primer sequence. This PCR pool is then amplified with exponential PCR, using primers that specifically amplify only the linear PCR molecules. An index primer is added to the amplicons along with a primer specific for each sample. This exponential PCR pool is then sequenced redundantly (20x coverage). The resulting sequences are separated by sample, using the index sequence. The amplicon sequences within each sample are separated by the random barcodes. The large number of reads for each barcode helps to create an error-corrected consensus sequence for the initial template molecule.
"""

script_info['script_usage'] = """Example: %prog -i fwd.fq,rev.fq -m Mapping_file.txt -o output_dir --barcode_type=7
"""

script_info['output_description'] = """ The %prog generates:\
A fasta file called seqs.fna which contains error corrected consensus sequence for the template DNA\
"""
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
                default=1.5),
    make_option('--min_consensus', type='float',
                help='threshold for consensus score'
                'the minimum score allowable at any position in sequence'
                'where the score is calulated as:'
                'occurence of base in consensus sequence/ total sequences'
                '[default: %default]',
                default=0.66),
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
	output_dir = opts.output_dir
	mapping_fp = opts.mapping_fp
	sequence_read_fps = opts.sequence_read_fps
	min_consensus = opts.min_consensus
	max_cluster_ratio = opts.threshold_for_cluster_ratio
	random_bc_lookup = read_input_file(sequence_read_fps, mapping_fp,\
                                       output_dir, barcode_type,\
                                       max_barcode_errors)
	consensus_outfile = open(os.path.join(output_dir, "seqs.fna"), "w")
	for sample_id in random_bc_lookup:
		for random_bc in random_bc_lookup[sample_id]:
			seqs, counts = random_bc_lookup[sample_id][random_bc]
			if cluster_ratio[sample_id][random_bc] < max_cluster_ratio:
				consensus_seq = select_majority_sequence(seqs, counts)
			else:
				consensus_seq = get_consensus(seqs, counts, min_consensus)  
			consensus_outfile.write(">"+ sample_id + random_bc 
                                        + "\n" + consensus_seq+"\n")
	
if __name__ == "__main__":
    main()
