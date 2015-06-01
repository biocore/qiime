#!/usr/bin/env python
# File created Sept 29, 2010
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from os.path import basename

from skbio.parse.sequences import parse_fasta

from qiime.parse import parse_qual_score


def parse_fasta_file(fasta_lines):
    """ Parses fasta file, generates dict of label:seq, list of seqs order


    fasta_lines: list of lines from fasta file.
    """

    fasta_seqs = {}
    seq_order = []

    for label, seq in parse_fasta(fasta_lines):
        fasta_seqs[label.split()[0].strip()] = seq
        seq_order.append(label)

    return fasta_seqs, seq_order


def verify_equivalency(fasta_seqs,
                       qual_scores):
    """ Tests for equivalent labels, base positions between fasta and qual file

    fasta_seqs:  dict of label:seq from fasta file
    qual_scores: dict of label: qual scores
    """

    if len(fasta_seqs) != len(qual_scores):
        raise ValueError('Number of sequences not equal in input fasta ' +
                         'and qual file.')

    qual_scores_labels = set(qual_scores.keys())

    for label in fasta_seqs.keys():

        # Should have equivalent labels
        if label not in qual_scores_labels:
            raise ValueError('Fasta label %s not found in quality score ' %
                             label + 'file.')

        # should have equivalent lengths
        if len(fasta_seqs[label]) != len(qual_scores[label]):
            raise ValueError('Sequence %s does not have equivalent ' %
                             label + 'base positions between fasta and quality score file.')


def truncate_seqs(fasta_seqs,
                  qual_scores,
                  base_pos):
    """ Truncates sequences to base position specified with base_pos

    fasta_seqs: dict of seq label: seq string
    qual_scores: dict of seq label: numpy array of int scores
    base_pos: index in sequence to truncate at
    """

    trunc_fasta_seqs = {}
    trunc_qual_scores = {}

    for seq in fasta_seqs:
        trunc_fasta_seqs[seq] = fasta_seqs[seq][:base_pos]
        trunc_qual_scores[seq] = qual_scores[seq][:base_pos]

    return trunc_fasta_seqs, trunc_qual_scores


def get_output_filepaths(output_dir,
                         fasta_fp,
                         qual_fp):
    """ Returns output filepaths for filtered fasta and quality files

    output_dir: output directory
    fasta_fp: input fasta filepath
    qual_fp: input quality scores filepath
    """

    if not output_dir.endswith('/'):
        output_dir += '/'

    fasta_out_fp = output_dir + basename(fasta_fp).split('.')[0] +\
        "_filtered.fasta"

    qual_out_fp = output_dir + basename(qual_fp).split('.')[0] +\
        "_filtered.qual"

    return fasta_out_fp, qual_out_fp


def write_trunc_fasta(trunc_fasta_seqs,
                      fasta_out_fp,
                      seq_order):
    """ Writes truncated fasta seqs in order specified with seq_order

    trunc_fasta_seqs: dict of fasta label: truncated sequence string
    fasta_out_fp: output filepath to write to
    seq_order: list of fasta labels in the order of the original input fasta
    """

    fasta_out = open(fasta_out_fp, "w")

    for label in seq_order:
        trunc_label = label.split()[0].strip()
        fasta_out.write(">%s\n%s\n" % (label, trunc_fasta_seqs[trunc_label]))


def write_trunc_qual(trunc_qual_scores,
                     qual_out_fp,
                     seq_order):
    """ Writes truncated quality score files out in proper format

    trunc_qual_scores: dict of seq label: numpy array of scores as ints
    qual_out_fp: output filepath to write truncated quality scores to
    seq_order: List of full fasta labels to write to output filepath and
     maintain the same order as input quality file.
    """

    qual_line_size = 60

    qual_out = open(qual_out_fp, "w")

    for label in seq_order:
        trunc_label = label.split()[0].strip()
        current_trunc_qual_scores = trunc_qual_scores[trunc_label]
        qual_out.write(">%s\n" % label)
        current_qual_scores_lines = []
        # Quality score format is a string of 60 base calls, followed by a
        # newline, until the last N bases are written
        for slice in range(0, len(trunc_qual_scores[trunc_label]),
                           qual_line_size):
            # current_segment = map(str,
            # current_trunc_qual_scores[slice:slice + qual_line_size])
            current_segment = current_trunc_qual_scores[
                slice:slice +
                qual_line_size]
            current_qual_scores_lines.append(" ".join(current_segment))

        qual_out.write('\n'.join(current_qual_scores_lines))
        qual_out.write('\n')


def truncate_fasta_qual(fasta_fp,
                        qual_fp,
                        output_dir,
                        base_pos):
    """ Main program function for generating quality score histogram

    fasta_fp: fasta filepath
    qual_fp: quality score filepath
    output_dir: output directory
    base_pos: Nucleotide position to truncate the fasta and quality score at.
    """

    qual_lines = open(qual_fp, "U")
    fasta_lines = open(fasta_fp, "U")

    qual_scores = parse_qual_score(qual_lines, value_cast_f=str)

    # Get dict of fasta label:seq, and the sequence order (so output can
    # be in the same order as the input sequences.
    fasta_seqs, seq_order = parse_fasta_file(fasta_lines)

    # Make sure the quality scores and fasta sequences have corresponding
    # labels and base numbers
    verify_equivalency(fasta_seqs, qual_scores)

    # Truncate seqs to base_pos index
    trunc_fasta_seqs, trunc_qual_scores = truncate_seqs(fasta_seqs,
                                                        qual_scores, base_pos)

    # Get output filepaths
    fasta_out_fp, qual_out_fp = get_output_filepaths(output_dir, fasta_fp,
                                                     qual_fp)

    # Write truncated sequences out
    write_trunc_fasta(trunc_fasta_seqs, fasta_out_fp, seq_order)

    write_trunc_qual(trunc_qual_scores, qual_out_fp, seq_order)
