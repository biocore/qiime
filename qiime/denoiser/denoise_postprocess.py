 #!/usr/bin/env python

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from optparse import OptionParser
from itertools import imap
from re import compile, search

from skbio.sequence import BiologicalSequence
from skbio.parse.sequences import parse_fasta

from qiime.denoiser.utils import read_denoiser_mapping, sort_ids


def extract_read_to_sample_mapping(labels):
    """Extract a mapping from reads to sample_ids from split_libraries.

    labels: iterable of header strings

      A fasta header from split_libraries looks like this
    >S160_1 E86FECS01DW5V4 orig_bc=CAGTACGATCTT new_bc=CAGTACGATCTT bc_diffs=0
    """
    sample_id_mapping = {}

    re = compile(r'(\S+) (\S+)')
    for label in labels:
        tmatch = search(re, label)
        sample_id = tmatch.group(1)
        flowgram_id = tmatch.group(2)
        sample_id_mapping[flowgram_id] = sample_id

    return sample_id_mapping


def post_process(fasta_fp, mapping_fp, denoised_seqs_fp,
                 otu_picker_otu_map_fp, out_dir):
    """redirect function for compatibilty with Qiime"""

    combine_mappings(open(fasta_fp), open(mapping_fp),
                     open(denoised_seqs_fp),
                     open(otu_picker_otu_map_fp), out_dir)


def combine_mappings(fasta_fh, mapping_fh, denoised_seqs_fh,
                     otu_picker_otu_map_fh, out_dir):
    """Combine denoiser and OTU picker mapping file, replace flowgram IDs.

    fasta_fh: a fasta file with labels as produced by Qiime's split_libraries.py
             used to replace flowgram id with the unique se_sample_id

    mapping_fh: The cluster mapping from the denoiser.py

    denoised_seqs_fh: the Fasta output files from denoiser.py

    otu_picker_map_fh: cluster map from otu picker on denoised_seqs_fh

    out_dir: output directory
    """

     # read in mapping from split_library file
    labels = imap(lambda a_b: a_b[0], parse_fasta(fasta_fh))
    # mapping from seq_id to sample_id
    sample_id_mapping = extract_read_to_sample_mapping(labels)

    denoiser_mapping = read_denoiser_mapping(mapping_fh)
    # read in cd_hit otu map
    # and write out combined otu_picker+denoiser map
    otu_fh = open(out_dir + "/denoised_otu_map.txt", "w")
    for otu_line in otu_picker_otu_map_fh:
        otu_split = otu_line.split()

        otu = otu_split[0]
        ids = otu_split[1:]

        get_sample_id = sample_id_mapping.get
        # concat lists
        # make sure the biggest one is first for pick_repr
        all_ids = sort_ids(ids, denoiser_mapping)
        all_ids.extend(sum([denoiser_mapping[id] for id in ids], []))
        try:
            otu_fh.write("%s\t" % otu +
                         "\t".join(map(get_sample_id, all_ids)) + "\n")
        except TypeError:
            # get returns Null if denoiser_mapping id not present in
            # sample_id_mapping
            print "Found id in denoiser output, which was not found in split_libraries " +\
                "output FASTA file. Wrong file?"
            exit()

    fasta_out_fh = open(out_dir + "/denoised_all.fasta", "w")
    for label, seq in parse_fasta(denoised_seqs_fh):
        id = label.split()[0]
        newlabel = "%s %s" % (sample_id_mapping[id], id)
        fasta_out_fh.write(BiologicalSequence(seq, id=newlabel).to_fasta())
