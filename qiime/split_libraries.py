#!/usr/bin/env python
# file split_libraries.py

"""Performs preprocessing steps for barcoded library analysis, e.g. 454.

Specifically, does the quality-filtering step (using several criteria) and
renames each read with the appropriate library id.

This module reads the linker+primer sequence from the input mapping file, and
associates these with the barcodes from the mapping file.  If a barcode is
read that does not correspond to any in the mapping file, this module checks
against all possible primers from the mapping file.  A rare situation could
arise if a barcode does not match any from the mapping file (either because
of sequencing errors or because the mapping file is incomplete) and variations
of the same primer are used for sequencing (e.g., a truncated form of the same
primer), where it is impossible to distinguish what primer was actually used
to amplify a given sequence.  In these cases, portions of the given sequence
are sliced out in ascending order of primer sizes and compared to all possible
primers from the mapping file.  The first matching primer is considered a hit
for the purposes of the determining where a primer ends and the actual
sequence data begins.  Because of this, one should be careful about using
truncated forms of the same primer with this module.  The -c option can be
used to disable attempts at barcode correction, and sequences that do not
have a barcode that matches the mapping file will not be recorded.
"""

__author__ = "Rob Knight and Micah Hamady"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Micah Hamady", "Greg Caporaso", "Kyle Bittinger",
               "Jesse Stombaugh", "William Walters", "Jens Reeder",
               "Emily TerAvest", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "rob@spot.colorado.edu, william.a.walters@colorado.edu"

import re
from gzip import GzipFile
from os import mkdir, stat
from collections import defaultdict
from string import upper

from numpy import array, mean, arange, histogram
from numpy import __version__ as numpy_version
import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from skbio.parse.sequences import parse_fasta
from cogent import DNA as DNA_cogent, LoadSeqs

from cogent.align.align import make_dna_scoring_dict, local_pairwise
from skbio.util import remove_files
from skbio.sequence import DNASequence

from qiime.check_id_map import process_id_map
from qiime.barcode import correct_barcode
from qiime.hamming import decode_barcode_8
from qiime.golay import decode as decode_golay_12
from qiime.format import format_histograms
from qiime.parse import QiimeParseError, parse_qual_scores
from qiime.util import create_dir, median_absolute_deviation

# Including new=True in the histogram() call is necessary to
# get the correct result in versions prior to NumPy 1.2.0,
# but the new keyword will be removed in NumPy 1.4. In NumPy 1.2.0
# or later, new=True raises a Warning regarding
# deprecation of new. One solution to this is to check the
# numpy version, and if it's less than 1.2.0, overwrite histogram
# with new=True as the default. This avoids a deprecation warning message
# in versions 1.2.0 through 1.3.*, and a try/except to handle errors from
# versions 1.4.0 or later.

numpy_version = re.split("[^\d]", numpy_version)
numpy_version = tuple([int(i) for i in numpy_version if i.isdigit()])
if numpy_version < (1, 3, 0):
    numpy_histogram = histogram

    def histogram(a, bins=10, range=None, normed=False, weights=None):
        return numpy_histogram(a, bins=bins, range=range,
                               normed=normed, weights=weights, new=True)

# Supported barcode types - need to adapt these functions to ignore list
# of valid barcodes that the generic decoder requires
BARCODE_TYPES = {
    "golay_12": (12, lambda bc, bcodes: decode_golay_12(bc)),
    "hamming_8": (8, lambda bc, bcodes: decode_barcode_8(bc)),
    # The decode function for variable length barcode does nothing -
    # it's just provided to comply with the interface of the other
    # barcode types. The returned barcode is always the same as the
    # one passed in, and the number of mismatches is always 0. The
    # length is None, corresponding to variable length.
    "variable_length": (None, lambda bc, bcodes: (bc, 0))}


def get_infile(filename):
    """Returns filehandle, allowing gzip input."""
    if filename.endswith(".gz"):
        fin = GzipFile(filename, "rb")
    else:
        fin = open(filename, "U")
    return fin


def count_mismatches(seq1, seq2, max_mm):
    """Counts mismatches, primer should be <= length of the seq.
    """
    mm = 0
    for i in range(min(len(seq2), len(seq1))):
        if seq1[i] != seq2[i]:
            mm += 1
            if mm > max_mm:
                return mm
    return mm


def ok_mm_primer(primer_seq, all_primer_seqs, primer_mm):
    """Check if primer_seq matches any primer within max mismatches.

    TODO: if we start using highly degenerate primers, should refactor using
    faster algorithm.
    """
    for curr_pat in all_primer_seqs:
        if count_mismatches(primer_seq, curr_pat, primer_mm) <= primer_mm:
            return True
    return False


def MatchScorerAmbigs(match, mismatch, matches=None):
    """ Alternative scorer factory for sw_align which allows match to ambiguous chars

    It allows for matching to ambiguous characters which is useful for
     primer/sequence matching. Not sure what should happen with gaps, but they
     shouldn't be passed to this function anyway. Currently a gap will only match
     a gap.

    match and mismatch should both be numbers. Typically, match should be
    positive and mismatch should be negative.

    Resulting function has signature f(x,y) -> number.

    Code original from Greg Caporaso
    """

    matches = matches or \
        {'A': {'A': None}, 'G': {'G': None}, 'C': {'C': None},
         'T': {'T': None}, '-': {'-': None}}
    for ambig, chars in DNASequence.iupac_degeneracies().iteritems():
        try:
            matches[ambig].update({}.fromkeys(chars))
        except KeyError:
            matches[ambig] = {}.fromkeys(chars)

        for char in chars:
            try:
                matches[char].update({ambig: None})
            except KeyError:
                matches[char] = {ambig: None}

    def scorer(x, y):
        # need a better way to disallow unknown characters (could
        # try/except for a KeyError on the next step, but that would only
        # test one of the characters)
        if x not in matches or y not in matches:
            raise ValueError("Unknown character: %s or %s" % (x, y))
        if y in matches[x]:
            return match
        else:
            return mismatch
    return scorer

# The scoring function which can be passed to
# cogent.alignment.algorithms.sw_align
equality_scorer_ambigs = MatchScorerAmbigs(1, -1)
expanded_equality_scorer_ambigs = MatchScorerAmbigs(1, -1,
                                                    matches=
                                                    {'A': {'A': None, 'G': None},
                                                     'G':
                                                     {'G': None,
                                                      'A': None,
                                                      'T': None,
                                                      'C': None,
                                                      },
                                                        'C':
                                                     {'C': None,
                                                      'G': None},
                                                        'T':
                                                     {'T': None,
                                                      'G': None},
                                                        '-': {'-': None}})


def pair_hmm_align_unaligned_seqs(seqs, moltype=DNA_cogent, params={}):
    """
        Checks parameters for pairwise alignment, returns alignment.

        Code from Greg Caporaso.
    """

    seqs = LoadSeqs(data=seqs, moltype=moltype, aligned=False)
    try:
        s1, s2 = seqs.values()
    except ValueError:
        raise ValueError(
            "Pairwise aligning of seqs requires exactly two seqs.")

    try:
        gap_open = params['gap_open']
    except KeyError:
        gap_open = 5
    try:
        gap_extend = params['gap_extend']
    except KeyError:
        gap_extend = 2
    try:
        score_matrix = params['score_matrix']
    except KeyError:
        score_matrix = make_dna_scoring_dict(
            match=1, transition=-1, transversion=-1)

    return local_pairwise(s1, s2, score_matrix, gap_open, gap_extend)


def local_align_primer_seq(primer, sequence, sw_scorer=equality_scorer_ambigs):
    """Perform local alignment of primer and sequence

        primer: Input primer sequence
        sequence: target sequence to test primer against

        Returns the number of mismatches,
         and the start position in sequence of the hit.

        Modified from code written by Greg Caporaso.
    """

    query_primer = primer

    query_sequence = str(sequence)

    # Get alignment object from primer, target sequence
    alignment = pair_hmm_align_unaligned_seqs([query_primer, query_sequence])

    # Extract sequence of primer, target site, may have gaps if insertions
    # or deletions have occurred.
    primer_hit = str(alignment.Seqs[0])
    target_hit = str(alignment.Seqs[1])

    # Count insertions and deletions
    insertions = primer_hit.count('-')
    deletions = target_hit.count('-')

    mismatches = 0
    for i in range(len(target_hit)):
        # using the scoring function to check for
        # matches, but might want to just access the dict
        if sw_scorer(target_hit[i], primer_hit[i]) == -1 and \
                target_hit[i] != '-' and primer_hit[i] != '-':
            mismatches += 1
    try:
        hit_start = query_sequence.index(target_hit.replace('-', ''))
    except ValueError:
        raise ValueError(
            'substring not found, query string %s, target_hit %s' %
            (query_sequence, target_hit))

    # sum total mismatches
    mismatch_count = insertions + deletions + mismatches

    return mismatch_count, hit_start


def expand_degeneracies(raw_primers):
    """Returns all non-degenerate versions of a given primer sequence.

    Order is not guaranteed!
    """

    expanded_primers = []

    for raw_primer in raw_primers:
        primer_seq = DNASequence(raw_primer.strip())

        for expanded_primer in primer_seq.nondegenerates():
            expanded_primers.append(str(expanded_primer))

    return expanded_primers


def check_map(infile, disable_primer_check, barcode_type="golay_12",
              added_demultiplex_field=None, has_barcodes=True):
    """Check mapping file and extract list of valid barcodes, primers """

    if barcode_type == "variable_length":
        var_len_barcodes = True
    else:
        var_len_barcodes = False

    if barcode_type == "0":
        has_barcodes = False

    # hds, id_map, dsp, run_description, errors, warnings
    hds, mapping_data, run_description, errors, warnings = \
        process_id_map(infile, has_barcodes=has_barcodes,
                       disable_primer_check=disable_primer_check,
                       added_demultiplex_field=added_demultiplex_field,
                       variable_len_barcodes=var_len_barcodes)

    if errors:
        raise ValueError('Errors were found with mapping file, ' +
                         'please run validate_mapping_file.py to ' +
                         'identify problems.')

    id_map = {}

    for curr_data in mapping_data:
        id_map[curr_data[0]] = {}

    for header in range(len(hds)):
        for curr_data in mapping_data:
            id_map[curr_data[0]][hds[header]] = curr_data[header]

    barcode_to_sample_id = {}

    primer_seqs_lens = {}
    all_primers = {}

    for sample_id, sample in id_map.items():
        if added_demultiplex_field:
            barcode_to_sample_id[sample['BarcodeSequence'].upper() + "," +
                                 sample[added_demultiplex_field]] = sample_id
        else:
            barcode_to_sample_id[sample['BarcodeSequence'].upper()] = sample_id
        if not disable_primer_check:
            raw_primers = sample['LinkerPrimerSequence'].upper().split(',')

            if len(raw_primers[0].strip()) == 0:
                raise ValueError('No primers detected, please use the ' +
                                 '-p parameter to disable primer detection.')
            expanded_primers = expand_degeneracies(raw_primers)
            curr_bc_primers = {}
            for primer in expanded_primers:
                curr_bc_primers[primer] = len(primer)
                all_primers[primer] = len(primer)
            primer_seqs_lens[sample['BarcodeSequence']] = curr_bc_primers

    return hds, id_map, barcode_to_sample_id, warnings, errors, \
        primer_seqs_lens, all_primers


def fasta_ids(fasta_files, verbose=False):
    """ Returns list of ids in FASTA files """
    all_ids = set([])
    for fasta_in in fasta_files:
        for label, seq in parse_fasta(fasta_in):
            rid = label.split()[0]
            if rid in all_ids:
                raise ValueError(
                    "Duplicate ID found in FASTA/qual file: %s" %
                    label)
            all_ids.add(rid)
    return all_ids


def count_ambig(curr_seq, valid_chars='ATCG'):
    """Counts non-standard characters in seq"""
    up_seq = curr_seq.upper()
    total = 0
    for vchar in valid_chars:
        total += up_seq.count(vchar)
    return len(curr_seq) - total


def split_seq(curr_seq, barcode_len, primer_seq_len):
    """Split sequence into parts barcode, primer and remainder"""
    curr_barcode = curr_seq[0:barcode_len]
    rest_of_seq = curr_seq[barcode_len:]
    primer_seq = rest_of_seq[0:primer_seq_len]
    rest_of_seq = rest_of_seq[primer_seq_len:]
    return curr_barcode, primer_seq, rest_of_seq


def get_barcode(curr_seq, barcode_len):
    """ Split sequence into barcode and remaining sequence

    Linker and primer part of remaining sequence, as one must first
    read the barcode to find the associated primer from the mapping file"""
    raw_barcode = curr_seq[0:barcode_len]
    raw_seq = curr_seq[barcode_len:]
    return raw_barcode, raw_seq


def primer_exceeds_mismatches(primer_seq, all_primer_seqs, max_primer_mm):
    """Returns True if primer exceeds allowed mismatches"""
    if primer_seq not in all_primer_seqs:
        if not ok_mm_primer(primer_seq, all_primer_seqs, max_primer_mm):
            return True
    return False


def seq_exceeds_homopolymers(curr_seq, max_len=6):
    """Returns False if primer contains any homopolymer > allowed length"""
    for base in 'ATGC':
        curr = base * (max_len + 1)
        if curr in curr_seq:
            return True
    return False


def check_barcode(curr_barcode, barcode_type, valid_map,
                  attempt_correction=True, added_demultiplex_field=None,
                  curr_id=None):
    """Return whether barcode is valid, and attempt correction."""

    corrected_bc = False

    if added_demultiplex_field:

        added_demultiplex_lens =\
            set([len(bc.split(',')[1]) for bc in valid_map])

        # using set() will put in order of smallest to largest and removes
        # redundant lengths, converting to list to sort from largest to
        # smallest
        added_demultiplex_lens =\
            [length for length in added_demultiplex_lens][::-1]

        # Handle specific case of run_prefix
        # Need to slice out size(s) of label that matches run prefix size(s)
        if added_demultiplex_field.upper() == "RUN_PREFIX":
            added_demultiplex =\
                [curr_id.split()[0][0:added_demultiplex_len] for
                 added_demultiplex_len in added_demultiplex_lens]

        else:
            for label_item in curr_id.split():
                if label_item.startswith(added_demultiplex_field):
                    added_demultiplex = [label_item.split('=')[1]]

        all_bcs = [bc.split(',')[0] for bc in valid_map]
        all_added_demultiplex = [bc.split(',')[1] for bc in valid_map]

        for curr_added_demultiplex in added_demultiplex:

            bc_and_demultiplex = curr_barcode + "," + curr_added_demultiplex

            if bc_and_demultiplex in valid_map:
                return False, bc_and_demultiplex, corrected_bc
            elif attempt_correction == False:
                return True, curr_barcode, corrected_bc
    else:
        if curr_barcode in valid_map:
            return False, curr_barcode, corrected_bc
        elif attempt_correction == False:
            return True, curr_barcode, corrected_bc

    if barcode_type in BARCODE_TYPES:
        expect_len, curr_bc_fun = BARCODE_TYPES[barcode_type]
        barcode, num_errors = curr_bc_fun(curr_barcode, valid_map)
        corrected_bc = True

        if added_demultiplex_field:

            for curr_added_demultiplex in added_demultiplex:
                bc_and_demultiplex = barcode + "," + curr_added_demultiplex
                if bc_and_demultiplex in valid_map:
                    return num_errors, bc_and_demultiplex, corrected_bc

        else:
            return num_errors, barcode, corrected_bc
    else:
        try:
            expect_len, curr_bc_fun = int(barcode_type), correct_barcode
            barcode, num_errors = curr_bc_fun(curr_barcode, valid_map)
            corrected_bc = True

            if added_demultiplex_field:

                for curr_added_demultiplex in added_demultiplex:
                    bc_and_demultiplex = barcode + "," + curr_added_demultiplex
                    if bc_and_demultiplex in valid_map:
                        return num_errors, bc_and_demultiplex, corrected_bc

        except ValueError:
            raise ValueError("Unsupported barcode type: %s" % barcode_type)
        return num_errors, barcode, corrected_bc


def make_histograms(raw_lengths, pre_lengths, post_lengths, binwidth=10):
    """Makes histogram data for pre and post lengths"""
    if post_lengths:
        min_len = min([min(post_lengths), min(raw_lengths)])
    else:
        min_len = min(raw_lengths)
    max_len = max(raw_lengths)
    floor = (min_len / binwidth) * binwidth
    ceil = ((max_len / binwidth) + 2) * binwidth
    bins = arange(floor, ceil, binwidth)
    raw_hist = histogram(raw_lengths, bins)[0]
    pre_hist = histogram(pre_lengths, bins)[0]
    post_hist, bin_edges = histogram(post_lengths, bins)
    return raw_hist, pre_hist, post_hist, bin_edges


class SeqQualBad(object):

    """Checks if a seq and qual score are bad, saving ids that are bad."""

    def __init__(self, name, f):
        """New SeqQualBad keeps track of failed ids."""
        self.FailedIds = []
        self.Name = name
        self.F = f

    def __call__(self, id_, seq, qual):
        """SeqQualBad called on id, seq and qual returns bool.

        Note: saves failed ids in self.FailedIds."""
        result = self.F(id_, seq, qual)
        if result:
            self.FailedIds.append(id_)
        return result

    def __str__(self):
        """SeqQualBad str returns tab-delimited output of counts."""
        return "%s\t%s" % (self.Name, len(self.FailedIds))


def qual_missing(id_, seq, qual):
    """returns True if qual is None"""
    return qual is None

QualMissing = SeqQualBad('Missing Qual Score', qual_missing)


def get_seq_lengths(seq_lengths, bc_counts):
    """Convenience wrapper for getting lengths of good and bad seqs"""
    all_seq_lengths = seq_lengths.values()
    all_seq_ids = set(seq_lengths.keys())
    bad_seq_ids = set(bc_counts[None]).union(set(bc_counts['#FAILED']))
    good_seq_ids = all_seq_ids - bad_seq_ids
    good_seq_lengths = map(seq_lengths.__getitem__, good_seq_ids)
    return all_seq_lengths, good_seq_lengths


def check_window_qual_scores(qual_scores, window=50, min_average=25):
    """Check that all windows have ave qual score > threshold."""

    # Code from Jens Reeder, added 1-13-2010
    l = len(qual_scores)

    window = min(window, l)
    if (window == 0):
        return True
    # initialize with sum of first window
    window_score = sum(qual_scores[:window])
    idx = 0
    while (window_score / float(window) >= min_average
           and idx < l - window):
            #'Move' window
        window_score += qual_scores[idx + window] - qual_scores[idx]
        idx += 1
    if (idx == l - window):
        # we processed all qual_scores, must be good
        # Return index for truncation purposes
        return True, idx
    else:
        return False, idx


def check_seqs(fasta_out, fasta_files, starting_ix, valid_map, qual_mappings,
               filters, barcode_len, keep_primer, keep_barcode, barcode_type,
               max_bc_errors, retain_unassigned_reads, attempt_bc_correction,
               primer_seqs_lens, all_primers, max_primer_mm, disable_primer_check,
               reverse_primers, rev_primers, qual_out, qual_score_window=0,
               discard_bad_windows=False, min_qual_score=25, min_seq_len=200,
               median_length_filtering=None, added_demultiplex_field=None,
               reverse_primer_mismatches=0, truncate_ambi_bases=False):
    """Checks fasta-format sequences and qual files for validity."""

    seq_lengths = {}

    # Record complete barcode + primer + sequence lengths
    raw_seq_lengths = {}
    # Record sequence lengths after all optional removal of components
    final_seq_lengths = {}

    bc_counts = defaultdict(list)
    curr_ix = starting_ix
    corr_ct = 0  # count of corrected barcodes

    # get the list of barcode lengths in reverse order
    barcode_length_order =\
        sorted(set([len(bc.split(',')[0]) for bc in valid_map]))
    barcode_length_order = barcode_length_order[::-1]

    primer_mismatch_count = 0
    all_primers_lens = sorted(set(all_primers.values()))

    reverse_primer_not_found = 0

    sliding_window_failed = 0
    trunc_ambi_base_counts = 0

    below_seq_min_after_trunc = 0
    below_seq_min_after_ambi_trunc = 0

    for fasta_in in fasta_files:
        for curr_id, curr_seq in parse_fasta(fasta_in):
            curr_rid = curr_id.split()[0]
            curr_seq = upper(curr_seq)

            curr_len = len(curr_seq)
            curr_qual = qual_mappings.get(curr_rid, None)

            # if qual_out:
            #    curr_qual_out_score = \
            #     "%2.2f" % float(float(sum(curr_qual))/float(len(curr_qual)))
            seq_lengths[curr_rid] = curr_len
            failed = False

            for f in filters:
                failed = failed or f(curr_rid, curr_seq, curr_qual)
            if failed:  # if we failed any of the checks, bail out here
                bc_counts['#FAILED'].append(curr_rid)
                continue

            if barcode_type == 'variable_length':
                # Reset the raw_barcode, raw_seq, and barcode_len -- if
                # we don't match a barcode from the mapping file, we want
                # these values to be None
                raw_barcode, raw_seq, barcode_len = (None, None, None)

                curr_valid_map =\
                    [curr_bc.split(',')[0] for curr_bc in valid_map]
                # Iterate through the barcode length from longest to shortest
                for l in barcode_length_order:
                    # extract the current length barcode from the sequence
                    bc, seq = get_barcode(curr_seq, l)
                    # check if the sliced sequence corresponds to a valid
                    # barcode, and if so set raw_barcode, raw_seq, and
                    # barcode_len for use in the next steps
                    if bc in curr_valid_map:
                        raw_barcode, raw_seq = bc, seq
                        barcode_len = len(raw_barcode)
                        break
                # if we haven't found a valid barcode, log this sequence as
                # failing to match a barcode, and move on to the next sequence
                if not raw_barcode:
                    bc_counts['#FAILED'].append(curr_rid)
                    continue

            else:
                # Get the current barcode to look up the associated primer(s)
                raw_barcode, raw_seq = get_barcode(curr_seq, barcode_len)

            if not disable_primer_check:
                try:
                    current_primers = primer_seqs_lens[raw_barcode]
                    # In this case, all values will be the same, i.e. the length
                    # of the given primer, or degenerate variations thereof.
                    primer_len = current_primers.values()[0]

                    if primer_exceeds_mismatches(raw_seq[:primer_len],
                                                 current_primers, max_primer_mm):
                        bc_counts['#FAILED'].append(curr_rid)
                        primer_mismatch_count += 1
                        continue
                except KeyError:
                    # If the barcode read does not match any of those in the
                    # mapping file, the situation becomes more complicated.  We do
                    # not know the length the sequence to slice out to compare to
                    # our primer sets, so, in ascending order of all the given
                    # primer lengths, a sequence will the sliced out and compared
                    # to the primer set.
                    current_primers = all_primers
                    found_match = False
                    for seq_slice_len in all_primers_lens:
                        if not(
                            primer_exceeds_mismatches(raw_seq[:seq_slice_len],
                                                      current_primers, max_primer_mm)):
                            primer_len = seq_slice_len
                            found_match = True
                            break
                    if not found_match:
                        bc_counts['#FAILED'].append(curr_rid)
                        primer_mismatch_count += 1
                        continue
                except IndexError:
                    # Try to raise meaningful error if problem reading primers
                    raise IndexError('Error reading primer sequences.  If ' +
                                     'primers were purposefully not included in the mapping ' +
                                     'file, disable usage with the -p option.')
            else:
                # Set primer length to zero if primers are disabled.
                primer_len = 0

            # split seqs
            cbc, cpr, cres = split_seq(curr_seq, barcode_len,
                                       primer_len)

            total_bc_primer_len = len(cbc) + len(cpr)

            # get current barcode
            try:
                bc_diffs, curr_bc, corrected_bc = \
                    check_barcode(cbc, barcode_type, valid_map.keys(),
                                  attempt_bc_correction, added_demultiplex_field, curr_id)
                if bc_diffs > max_bc_errors:
                    raise ValueError("Too many errors in barcode")
                corr_ct += bool(corrected_bc)
            except Exception as e:
                bc_counts[None].append(curr_rid)
                continue

            curr_samp_id = valid_map.get(curr_bc, 'Unassigned')

            new_id = "%s_%d" % (curr_samp_id, curr_ix)
            # check if writing out primer
            write_seq = cres

            if reverse_primers == "truncate_only":
                try:
                    rev_primer = rev_primers[curr_bc]
                    mm_tested = {}
                    for curr_rev_primer in rev_primer:
                        # Try to find lowest count of mismatches for all
                        # reverse primers
                        rev_primer_mm, rev_primer_index  = \
                            local_align_primer_seq(curr_rev_primer, cres)
                        mm_tested[rev_primer_mm] = rev_primer_index

                    rev_primer_mm = min(mm_tested.keys())
                    rev_primer_index = mm_tested[rev_primer_mm]
                    if rev_primer_mm <= reverse_primer_mismatches:
                        write_seq = write_seq[0:rev_primer_index]
                        if qual_out:
                            curr_qual = curr_qual[0:barcode_len +
                                                  primer_len + rev_primer_index]
                    else:
                        reverse_primer_not_found += 1
                except KeyError:
                    pass
            elif reverse_primers == "truncate_remove":
                try:
                    rev_primer = rev_primers[curr_bc]
                    mm_tested = {}
                    for curr_rev_primer in rev_primer:
                        # Try to find lowest count of mismatches for all
                        # reverse primers
                        rev_primer_mm, rev_primer_index  = \
                            local_align_primer_seq(curr_rev_primer, cres)
                        mm_tested[rev_primer_mm] = rev_primer_index

                    rev_primer_mm = min(mm_tested.keys())
                    rev_primer_index = mm_tested[rev_primer_mm]
                    if rev_primer_mm <= reverse_primer_mismatches:
                        write_seq = write_seq[0:rev_primer_index]
                        if qual_out:
                            curr_qual = curr_qual[0:barcode_len +
                                                  primer_len + rev_primer_index]
                    else:
                        reverse_primer_not_found += 1
                        write_seq = False
                except KeyError:
                    bc_counts['#FAILED'].append(curr_rid)
                    continue

            # Check for quality score windows, truncate or remove sequence
            # if poor window found.  Previously tested whole sequence-now
            # testing the post barcode/primer removed sequence only.
            if qual_score_window:
                passed_window_check, window_index =\
                    check_window_qual_scores(curr_qual, qual_score_window,
                                             min_qual_score)
                # Throw out entire sequence if discard option True
                if discard_bad_windows and not passed_window_check:
                    sliding_window_failed += 1
                    write_seq = False
                # Otherwise truncate to index of bad window
                elif not discard_bad_windows and not passed_window_check:
                    sliding_window_failed += 1
                    if write_seq:
                        write_seq = write_seq[0:window_index]
                        if qual_out:
                            curr_qual = curr_qual[0:barcode_len +
                                primer_len + window_index]
                        #Check for sequences that are too short after truncation
                        if len(write_seq) + total_bc_primer_len < min_seq_len:
                            write_seq = False
                            below_seq_min_after_trunc += 1

            if truncate_ambi_bases and write_seq:
                write_seq_ambi_ix = True
                # Skip if no "N" characters detected.
                try:
                    ambi_ix = write_seq.index("N")
                    write_seq = write_seq[0:ambi_ix]
                except ValueError:
                    write_seq_ambi_ix = False
                    pass
                if write_seq_ambi_ix:
                    # Discard if too short after truncation
                    if len(write_seq) + total_bc_primer_len < min_seq_len:
                        write_seq = False
                        below_seq_min_after_ambi_trunc += 1
                    else:
                        trunc_ambi_base_counts += 1
                        if qual_out:
                            curr_qual = curr_qual[0:barcode_len +
                                                  primer_len + ambi_ix]

            # Slice out regions of quality scores that correspond to the
            # written sequence, i.e., remove the barcodes/primers and reverse
            # primers if option is enabled.
            if qual_out:
                qual_barcode, qual_primer, qual_scores_out = \
                    split_seq(curr_qual, barcode_len, primer_len)
                # Convert to strings instead of numpy arrays, strip off
                # brackets
                qual_barcode = format_qual_output(qual_barcode)
                qual_primer = format_qual_output(qual_primer)
                qual_scores_out = format_qual_output(qual_scores_out)

            if not write_seq:
                bc_counts['#FAILED'].append(curr_rid)
                continue

            if keep_primer:
                write_seq = cpr + write_seq
                if qual_out:
                    qual_scores_out = qual_primer + qual_scores_out
            if keep_barcode:
                write_seq = cbc + write_seq
                if qual_out:
                    qual_scores_out = qual_barcode + qual_scores_out

            # Record number of seqs associated with particular barcode.
            bc_counts[curr_bc].append(curr_rid)

            if retain_unassigned_reads and curr_samp_id == "Unassigned":
                fasta_out.write(
                    ">%s %s orig_bc=%s new_bc=%s bc_diffs=%s\n%s\n" %
                    (new_id, curr_rid, cbc, curr_bc, int(bc_diffs), write_seq))
                if qual_out:
                    qual_out.write(
                        ">%s %s orig_bc=%s new_bc=%s bc_diffs=%s\n%s" %
                        (new_id, curr_rid, cbc, curr_bc, int(bc_diffs),
                         qual_scores_out))
            elif not retain_unassigned_reads and curr_samp_id == "Unassigned":
                bc_counts['#FAILED'].append(curr_rid)
            else:
                fasta_out.write(
                    ">%s %s orig_bc=%s new_bc=%s bc_diffs=%s\n%s\n" %
                    (new_id, curr_rid, cbc, curr_bc, int(bc_diffs), write_seq))
                if qual_out:
                    qual_out.write(
                        ">%s %s orig_bc=%s new_bc=%s bc_diffs=%s\n%s" %
                        (new_id, curr_rid, cbc, curr_bc, int(bc_diffs),
                         qual_scores_out))

            curr_len = len(write_seq)

            #seq_lengths[curr_rid] = curr_len

            curr_ix += 1

            # Record the raw and written seq length of everything passing
            # filters
            raw_seq_lengths[curr_rid] = len(curr_seq)
            final_seq_lengths[curr_id] = curr_len

    if median_length_filtering:
        # Read original fasta file output to get sequence lengths
        fasta_out.close()
        fasta_out = open(fasta_out.name, "U")

        # Record sequence lengths for median/mad calculation
        sequence_lens = []
        for label, seq in parse_fasta(fasta_out):
            sequence_lens.append(len(seq))

        '''# Create a temporary file to copy the contents of the fasta file, will
        # need to delete once operations complete.
        fasta_temp = open(fasta_out.name + "_tmp.fasta", "w")

        sequence_lens = []
        for label, seq in parse_fasta(fasta_lens):
            sequence_lens.append(len(seq))
            fasta_temp.write(">%s\n%s\n" % (label, seq))

        fasta_temp.close()
        fasta_temp = open(fasta_out.name + "_tmp.fasta", "U")

        fasta_lens.close()
        # Overwrite seqs.fna with length filtered data
        fasta_out = open(fasta_out.name, "w")'''

        med_abs_dev, med_length = median_absolute_deviation(sequence_lens)

        min_corrected_len = med_length - med_abs_dev *\
            float(median_length_filtering)
        max_corrected_len = med_length + med_abs_dev *\
            float(median_length_filtering)
        seqs_discarded_median = 0

        fasta_out.seek(0)

        final_written_lens = []

        # Create final seqs.fna
        final_fasta_out = open(fasta_out.name.replace('.tmp', ''), "w")

        for label, seq in parse_fasta(fasta_out):
            curr_len = len(seq)
            if curr_len < min_corrected_len or curr_len > max_corrected_len:
                seqs_discarded_median += 1
            else:
                final_fasta_out.write(">%s\n%s\n" % (label, seq))
                final_written_lens.append(len(seq))

        final_fasta_out.close()
        fasta_out.close()
        remove_files([fasta_out.name])

    else:
        min_corrected_len = 0
        max_corrected_len = 0
        seqs_discarded_median = 0
        final_written_lens = 0

        # Copy tmp seqs file to final seqs.fna file
        fasta_out.close()
        fasta_out = open(fasta_out.name, "U")

        # Create final seqs.fna
        final_fasta_out = open(fasta_out.name.replace('.tmp', ''), "w")

        for label, seq in parse_fasta(fasta_out):
            final_fasta_out.write(">%s\n%s\n" % (label, seq))

        final_fasta_out.close()
        fasta_out.close()
        remove_files([fasta_out.name])

    median_results = (median_length_filtering, min_corrected_len,
                      max_corrected_len, seqs_discarded_median, final_written_lens)

    raw_seq_lengths = raw_seq_lengths.values()
    final_seq_lengths = final_seq_lengths.values()

    log_out = format_log(bc_counts, corr_ct, valid_map, seq_lengths, filters,
                         retain_unassigned_reads, attempt_bc_correction, primer_mismatch_count,
                         max_primer_mm, reverse_primers, reverse_primer_not_found,
                         sliding_window_failed, below_seq_min_after_trunc, qual_score_window,
                         discard_bad_windows, min_seq_len, raw_seq_lengths,
                         final_seq_lengths, median_results, truncate_ambi_bases,
                         below_seq_min_after_ambi_trunc, )

    #all_seq_lengths, good_seq_lengths = get_seq_lengths(seq_lengths, bc_counts)

    return log_out, seq_lengths.values(), raw_seq_lengths, final_seq_lengths


def format_qual_output(qual_array):
    """ Converts to string from numpy arrays, removes brackets """

    # Size of lines needed for proper quality score file format
    qual_line_size = 60

    qual_scores = ""

    for slice in range(0, len(qual_array), qual_line_size):
        current_segment = qual_array[slice:slice + qual_line_size]
        current_segment =\
            " ".join(str(score) for score in current_segment) + "\n"

        qual_scores += current_segment

    '''qual_array = str(qual_array)
    qual_array = qual_array.replace('[','')
    qual_array = qual_array.replace(']','') '''
    return qual_scores


def format_log(bc_counts, corr_ct, valid_map, seq_lengths, filters,
               retain_unassigned_reads, attempt_bc_correction,
               primer_mismatch_count, max_primer_mm,
               reverse_primers, reverse_primer_not_found, sliding_window_failed,
               below_seq_min_after_trunc, qual_score_window,
               discard_bad_windows, min_seq_len,
               raw_seq_lengths, final_seq_lengths, median_results=(None),
               truncate_ambi_bases=False, below_seq_min_after_ambi_trunc=0,
               ):
    """Makes log lines"""
    log_out = []
    all_seq_lengths, good_seq_lengths = get_seq_lengths(seq_lengths, bc_counts)
    log_out.append("Number raw input seqs\t%d\n" % len(seq_lengths))

    # append log data for median absolute deviation sequence length filtering
    # if was performed.
    if median_results[0]:
        if (not median_results[1] or not median_results[2] or
                not median_results[3]):
            log_out.append("No sequences written, so no median length data " +
                           "available.")
            actual_median_results = False
        else:
            log_out.append("Specified allowed number of median absolute " +
                           "deviations for sequence retention: %3.2f" % (float(median_results[0])))
            log_out.append("Sequences with lengths outside bounds of " +
                           "%d and %d: %d\n" %
                           (int(median_results[1]), int(median_results[2]),
                            int(median_results[3])))
            actual_median_results = True

    for f in filters:
        log_out.append(str(f))
    log_out.append('Num mismatches in primer exceeds limit of %s: %d\n' %
                   (max_primer_mm, primer_mismatch_count))
    if reverse_primers == "truncate_only":
        log_out.append('Number of sequences with identifiable barcode ' +
                       'but without identifiable reverse primer: ' +
                       '%d\n' % reverse_primer_not_found)
        log_out.append('-z truncate_only option enabled; sequences ' +
                       'without a discernible reverse primer as well as sequences with a ' +
                       'valid barcode not found in the mapping file may still be written.\n')
    if reverse_primers == "truncate_remove":
        log_out.append('Number of sequences with identifiable barcode ' +
                       'but without identifiable reverse primer: ' +
                       '%d\n' % reverse_primer_not_found)
        log_out.append('-z truncate_remove option enabled; sequences ' +
                       'without a discernible reverse primer as well as sequences with a ' +
                       'valid barcode not found in the mapping file will not be written.\n')
    if qual_score_window:
        log_out.append('Size of quality score window, in base pairs: %d' %
                       qual_score_window)
        log_out.append('Number of sequences where a low quality score ' +
                       'window was detected: %d' % sliding_window_failed)
        if discard_bad_windows:
            log_out.append('Sequences with a low quality score were not ' +
                           'written, -g option enabled.\n')
        else:
            log_out.append('Sequences with low quality score window were ' +
                           'truncated to the first base of the window.')
            log_out.append('Sequences discarded after truncation due to ' +
                           'sequence length below the minimum %d: %d\n' %
                           (min_seq_len, below_seq_min_after_trunc))
    if truncate_ambi_bases:
        log_out.append('Truncation at first ambiguous "N" character ' +
                       'enabled.\nSequences discarded after truncation due to sequence ' +
                       'length below the minimum %d: %d\n' %
                       (min_seq_len, below_seq_min_after_ambi_trunc))

    log_out.append("Sequence length details for all sequences passing " +
                   "quality filters:")
    if raw_seq_lengths:
        log_out.append("Raw len min/max/avg\t%.1f/%.1f/%.1f" %
                       (min(raw_seq_lengths), max(raw_seq_lengths), mean(raw_seq_lengths)))
    else:
        log_out.append("No sequences passed quality filters for writing.")

    if median_results[0] and actual_median_results:
        log_out.append("Wrote len min/max/avg\t%.1f/%.1f/%.1f" %
                       (min(median_results[4]), max(median_results[4]),
                        mean(median_results[4])))
    else:
        if final_seq_lengths:
            log_out.append("Wrote len min/max/avg\t%.1f/%.1f/%.1f" %
                           (min(final_seq_lengths), max(final_seq_lengths),
                            mean(final_seq_lengths)))

    # figure out which barcodes we got that didn't come from valid samples
    valid_bc_nomap = set(bc_counts) - set(valid_map) - set([None, '#FAILED'])
    valid_bc_nomap_counts = [(len(bc_counts[b]), b) for b in valid_bc_nomap]

    log_out.append("\nBarcodes corrected/not\t%d/%d" %
                   (corr_ct, len(bc_counts[None])))
    if attempt_bc_correction:
        log_out.append("Uncorrected barcodes will not be written to the " +
                       "output fasta file.\nCorrected barcodes will be written with " +
                       "the appropriate barcode category.\nCorrected but unassigned " +
                       "sequences will not be written unless --retain_unassigned_reads " +
                       "is enabled.\n")
    else:
        log_out.append("Barcode correction has been disabled via the -c " +
                       "option.\n")
    log_out.append("Total valid barcodes that are not in mapping file\t%d" %
                   len(valid_bc_nomap_counts))
    if valid_bc_nomap:
        log_out.append("Barcodes not in mapping file\tCount")
        for count, bc in reversed(sorted(valid_bc_nomap_counts)):
            log_out.append("%s\t%d" % (bc, count))
    if retain_unassigned_reads:
        log_out.append("Sequences associated with valid barcodes that are not" +
                       " in the mapping file will be written. --retain_unassigned_reads " +
                       "option enabled.")
    elif not attempt_bc_correction:
        log_out.append("Barcode correction has been disabled (-c option), " +
                       "no unassigned or invalid barcode sequences will be recorded.")
    else:
        log_out.append("Sequences associated with valid barcodes that are " +
                       "not in the mapping file will not be written.")

    log_out.append("\nBarcodes in mapping file")

    sample_cts = [(len(bc_counts[bc]), bc, sample_id) for bc, sample_id
                  in valid_map.items()]
    if sample_cts:
        filtered_sample_cts = [s[0] for s in sample_cts if s[0]]
        if filtered_sample_cts:
            log_out.append("Num Samples\t%d" % len(filtered_sample_cts))
            log_out.append("Sample ct min/max/mean: %d / %d / %.2f" % (
                min(filtered_sample_cts), max(filtered_sample_cts),
                mean(filtered_sample_cts)))
        log_out.append("Sample\tSequence Count\tBarcode")
        for count, bc, sample_id in reversed(sorted(sample_cts)):
            log_out.append("%s\t%d\t%s" % (sample_id, count, bc))

    if median_results[0]:
        corrected_written_seqs = len(good_seq_lengths) - median_results[3]
    else:
        corrected_written_seqs = len(good_seq_lengths)

    log_out.append("\nTotal number seqs written\t%d" % corrected_written_seqs)
    return log_out


def get_reverse_primers(id_map):
    """ Return a dictionary with barcodes and rev-complement of rev primers """

    rev_primers = {}
    for n in id_map.items():
        # Generate a dictionary with Barcode:reverse primer
        # Convert to reverse complement of the primer so its in the
        # proper orientation with the input fasta sequences
        rev_primers[n[1]['BarcodeSequence']] =\
            [str(DNASequence(curr_rev_primer).rc()) for curr_rev_primer in
             (n[1]['ReversePrimer']).split(',')]

    return rev_primers


def preprocess(fasta_files, qual_files, mapping_file,
               barcode_type="golay_12",
               min_seq_len=200, max_seq_len=1000, min_qual_score=25, starting_ix=1,
               keep_primer=False, max_ambig=0, max_primer_mm=0, trim_seq_len=False,
               dir_prefix='.', max_bc_errors=2, max_homopolymer=4,
               retain_unassigned_reads=False, keep_barcode=False,
               attempt_bc_correction=True, qual_score_window=0,
               disable_primer_check=False, reverse_primers='disable',
               reverse_primer_mismatches=0,
               record_qual_scores=False, discard_bad_windows=False,
               median_length_filtering=None, added_demultiplex_field=None,
               truncate_ambi_bases=False):
    """
    Preprocess barcoded libraries, e.g. from 454.

    Parameters:

    fasta_files: list of raw 454 fasta files, fasta format.

    qual_files: list of raw 454 qual file(s)

    mapping_file: mapping file with BarcodeSequence column containing valid
    barcodes used in the 454 run

    barcode_type: type of barcode, e.g. golay_12. Should appear in list of
    known barcode types.

    min_seq_len: minimum sequence length to allow.

    max_seq_len: maximum sequence length to allow.

    min_qual_score: minimum average qual score considered acceptaable.

    starting_ix: integer to start sample sequence numbering at.

    keep_primer: when True, will keep primer sequence, otherwise will strip it

    keep_barcode: when True, will keep barcode sequence, otherwise will strip it

    max_ambig: maximum number of ambiguous bases to allow in the read.

    max_primer_mm: maximum number of primer mismatches to allow.

    trim_seq_len: if True, calculates lengths after trimming.

    dir_prefix: prefix of directories to write files into.

    max_bc_errors: maximum number of barcode errors to allow in output seqs

    max_homopolymer: maximum number of a nucleotide that can be
    repeated in a given sequence.

    retain_unassigned_reads: If True (False default), will write seqs to the
    output .fna file that have a valid barcode (by Golay or Hamming standard)
    but are not included in the input mapping file, as "Unassigned."

    attempt_bc_correction: (default True) will attempt to find nearest valid
    barcode.  Can be disabled to improve performance.

    disable_primer_check: (default False) Disables testing for primers in the
    input mapping file and primer testing in the input sequence files.

    reverse_primers: (default 'disable') Enables removal of reverse primers and
    any subsequence sequence data from output reads.  Reverse primers have to
    be in 5'->3' format and in correct IUPAC codes in a column "ReversePrimer"
    in the input mapping file.  Run validate_mapping_file.py to make test primers in this
    column for valid formatting.  The primers read from this column will be
    reverse complemented and associated with the given barcode in the
    mapping file.  If set to 'truncate_only', sequences where primers are found
    will be truncated, sequences where the primer is not found will be written
    unchanged.  If set to 'truncate_remove', sequences where primers are found
    will be truncated, sequences where the primer is not found will not be
    written and counted in the log file as failing for this reason.  The
    mismatches allowed for a reverse primer match are the same as specified
    for the forward primer mismatches with the -M parameter (default 0).

    reverse_primer_mismatches: Number of reverse primer mismatches allowed.
    reverse_primers must be enabled for this to do anything.

    record_qual_scores:  (default False) Will record quality scores for all
    sequences that are written to the output seqs.fna file in a separate
    file (seqs_filtered.qual) containing the same sequence IDs and
    quality scores for all bases found in the seqs.fna file.

    discard_bad_windows: (default False) If True, will completely discard
    sequences that have a low quality window.  If False, sequences will be
    truncated to the first base of the bad window.

    median_length_filtering: (default None) If a value is specified, will
    disable all min and max length filtering, and instead will calculate the
    median length of all sequences output, and instead filter out sequences
    based upon whether or not they fall within the number of median absolute
    deviations given by this parameter.

    added_demultiplex_field: (default None) If enabled, will attempt to
    demultiplex by both a barcode and data that can be parsed from the fasta
    label/comment.  If 'run_prefix' is specified will pull the character string
    immediatly following the '>' at the beginning of the fasta label.  Any
    other field specified will be pulled from the comment.  Example: if 'plate'
    is specified, the string following 'plate=' will be used.  The mapping file
    must contain a column with a header matching the name specified by the -j
    option, and every combination of barcode + added demultiplex option must
    be unique.

    truncate_ambi_bases: (default False) If enabled, will truncate the
    sequence at the first "N" character.

    Result:
    in dir_prefix, writes the following files:
    id_map.xls: 2-column tab-delimited text format orig_id:new_id
    error_map.xls: 2-column tab-delimited text format orig_id:fail_reasons
    seqs.fasta: sequences with new ids lib_index in fasta format
    lengths.xls: histograms of unfiltered and filtered lengths, resolution 10 bp
    """

    if max_seq_len < 10:
        raise ValueError("Max sequence must be >= 10")
    if min_seq_len >= max_seq_len:
        raise ValueError("Min len cannot be >= max len")
    if min_qual_score < 0:
        raise ValueError("Min qual score must be > 0")
    if starting_ix < 1:
        raise ValueError("Starting index must be > 0.")
    if max_ambig < 0:
        raise ValueError("Max ambig chars must be >= 0.")
    if max_primer_mm < 0:
        raise ValueError("Max primer mismatches must be >= 0.")
    if reverse_primers not in ['disable', 'truncate_only', 'truncate_remove']:
        raise ValueError("reverse_primers parameter must be 'disable', " +
                         "truncate_only, or truncate_remove.")

    create_dir(dir_prefix, fail_on_exist=False)

#    try:
#        stat(dir_prefix)
#    except OSError:
#        mkdir(dir_prefix)

    """# Generate primer sequence patterns - changing to mapping file primers.
    all_primer_seqs, primer_seq_len = \
        get_primer_seqs(primer_seq_pats.split(',')) """

    # Check mapping file and get barcode mapping
    map_file = open(mapping_file, 'U')
    headers, id_map, valid_map, warnings, errors, \
        primer_seqs_lens, all_primers = check_map(
            map_file, disable_primer_check,
            barcode_type, added_demultiplex_field)

    if reverse_primers != 'disable':
        if 'ReversePrimer' not in headers:
            raise ValueError('To enable reverse primer check, there must ' +
                             'be a "ReversePrimer" column in the mapping file with a reverse ' +
                             'primer in each cell.')
        rev_primers = get_reverse_primers(id_map)
    else:
        rev_primers = False

    # *** Generate dictionary of {barcode: DNA(ReversePrimer).rc()}
    # First check for ReversePrimer in headers, raise error if not found
    # Implement local alignment for primer after barcode is determined.
    # Add option to flag seq with error for rev_primer not found
    # Check primer hit index, truncate sequence
    # unit tests.

    map_file.close()
    if errors:
        raise ValueError("Invalid mapping file. " +
                         "Validate with validate_mapping_file.py first: %s" % "\n".join(errors))

    # Find actual length of barcodes in the mapping file, also check for
    # variable lengths, in case of added_demultiplex, split on comma.
    barcode_length_check =\
        list(set([len(bc.split(',')[0]) for bc in valid_map]))

        # Check barcode type
    if barcode_type not in BARCODE_TYPES:
        try:
            barcode_len, barcode_fun = int(barcode_type), correct_barcode
        except ValueError:
            raise ValueError("Unsupported barcode type: %s" % barcode_type)
    else:
        barcode_len, barcode_fun = BARCODE_TYPES[barcode_type]

    # As people often do not specify a barcode that matches the lengths
    # of the barcodes used, a check on the actual barcode lengths needs to
    # be done, and an exception raised if they are variable length and not
    # specified as so.
    if barcode_type != "variable_length":
        # Raise error if variable length barcodes are present but not
        # specified
        if len(barcode_length_check) != 1:
            raise ValueError('Mapping file has variable length ' +
                             'barcodes.  If this is intended, specifiy variable lengths ' +
                             'with the -b variable_length option.')
        # Raise error if the specified barcode length doesn't match what
        # is present in the mapping file.
        if barcode_len != barcode_length_check[0]:
            raise ValueError('Barcode length detected in the mapping file, ' +
                             ' %d does not match specified barcode length, %d.  ' %
                             (barcode_length_check[0], barcode_len) + 'To specify a barcode ' +
                             'length use -b golay_12 or -b hamming_8 for 12 and 8 base pair ' +
                             'golay or hamming codes respectively, or -b # where # is the ' +
                             'length of the barcode used.  E.g. -b 4 for 4 base pair barcodes.')

    fasta_files = map(get_infile, fasta_files)
    qual_files = map(get_infile, qual_files)

    # Check fasta files valid format, no duplicate ids
    # and ids match between fasta and qual files
    all_fasta_ids = fasta_ids(fasta_files)
    all_qual_ids = fasta_ids(qual_files)
    if qual_files and (len(all_fasta_ids) != len(all_qual_ids)):
        f_ids = all_fasta_ids.difference(all_qual_ids)
        q_ids = all_qual_ids.difference(all_fasta_ids)
        raise ValueError(
            "Found %d ids in fasta file not in qual file, %d ids in qual file not in fasta" %
            (len(f_ids), len(q_ids)))

    for f in fasta_files:
        f.seek(0)
    if qual_files:
        for q in qual_files:
            q.seek(0)
        # Load quality scores
        qual_mappings = parse_qual_scores(qual_files)
        for q in qual_files:
            q.close()
    else:
        qual_mappings = {}

    # make filters
    filters = []
    # seq len filter depends on whether we're including the barcode, if
    # median_length_filtering turned on, no length filtering.
    if not median_length_filtering:
        if trim_seq_len:
            # This processing occurs before primer testing, will use largest
            # primer length to calculate lengths.  the dict all_primers has
            # keys of each primer with the length of said primer as the value
            if disable_primer_check:
                primer_seq_len = 0
            else:
                primer_seq_len = max(all_primers.values())

            if barcode_type == "variable_length":
                barcode_len = max(barcode_length_check)

            trim = barcode_len + primer_seq_len
            filters.append(SeqQualBad(
                'Length outside bounds of %s and %s' % (
                    min_seq_len,
                    max_seq_len),
                lambda id_, seq, qual:
                not (min_seq_len <= len(seq) - trim <= max_seq_len)))
        else:
            filters.append(SeqQualBad(
                'Length outside bounds of %s and %s' % (
                    min_seq_len,
                    max_seq_len),
                lambda id_, seq, qual: not (min_seq_len <= len(seq) <= max_seq_len)))

    if not truncate_ambi_bases:
        filters.append(SeqQualBad(
            'Num ambiguous bases exceeds limit of %s' % max_ambig,
            lambda id_, seq, qual: count_ambig(seq) > max_ambig))

    if qual_mappings:
        filters.append(QualMissing)
        filters.append(SeqQualBad(
            'Mean qual score below minimum of %s' % min_qual_score,
            lambda id_, seq, qual: mean(qual) < min_qual_score))
    """if qual_score_window:
            filters.append(SeqQualBad('Mean window qual score below '+\
                            'minimum of %s' % min_qual_score,
                                      lambda id_, seq, qual: \
             not check_window_qual_scores(qual, qual_score_window, \
             min_qual_score))) """

    # Changed this to check entire sequence after barcode-could cause issue
    # if barcode-linker-primer have long homopolymers though.
    filters.append(SeqQualBad(
        'Max homopolymer run exceeds limit of %s' % max_homopolymer,
        lambda id_, seq, qual: seq_exceeds_homopolymers(
            seq[barcode_len:], max_homopolymer)))

    # Check seqs and write out
    fasta_out = open(dir_prefix + '/' + 'seqs.fna.tmp', 'w+')
    if record_qual_scores:
        qual_out = open(dir_prefix + '/' + 'seqs_filtered.qual', 'w+')
    else:
        qual_out = False

    '''log_stats, pre_lens, post_lens = check_seqs(fasta_out, fasta_files,
        starting_ix, valid_map, qual_mappings, filters, barcode_len,
        primer_seq_len, keep_primer, keep_barcode, barcode_type, max_bc_errors,
        retain_unassigned_reads) '''
    log_stats, raw_lens, pre_lens, post_lens = check_seqs(fasta_out,
                                                          fasta_files, starting_ix, valid_map, qual_mappings, filters,
                                                          barcode_len, keep_primer, keep_barcode, barcode_type, max_bc_errors,
                                                          retain_unassigned_reads, attempt_bc_correction,
                                                          primer_seqs_lens, all_primers, max_primer_mm, disable_primer_check,
                                                          reverse_primers, rev_primers, qual_out, qual_score_window,
                                                          discard_bad_windows, min_qual_score, min_seq_len,
                                                          median_length_filtering, added_demultiplex_field,
                                                          reverse_primer_mismatches, truncate_ambi_bases)

    # Write log file
    log_file = open(dir_prefix + '/' + "split_library_log.txt", 'w+')
    log_file.write('\n'.join(log_stats))
    log_file.close()

    # Write sequence distros here
    histogram_file = open(dir_prefix + '/' + 'histograms.txt', 'w+')

    histogram_file.write(format_histograms
                         (*make_histograms(raw_lens, pre_lens, post_lens)))
    histogram_file.close()
