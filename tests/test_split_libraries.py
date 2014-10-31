#!/usr/bin/env python
# file test_split_libraries.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "William Walters", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"


from os import close
from os.path import exists
from StringIO import StringIO
from numpy import array
from shutil import rmtree
from tempfile import mkstemp, mkdtemp

from unittest import TestCase, main
from numpy. testing import assert_almost_equal
from skbio.util import remove_files

from qiime.split_libraries import (
    expand_degeneracies, get_infile, count_mismatches,
    ok_mm_primer, check_map, fasta_ids,
    count_ambig, split_seq, primer_exceeds_mismatches,
    check_barcode, make_histograms, SeqQualBad,
    seq_exceeds_homopolymers, check_window_qual_scores, check_seqs,
    local_align_primer_seq, preprocess)
from qiime.parse import parse_qual_score

class FakeOutFile(object):

    def __init__(self):
        self.data = ""

    def write(self, s):
        self.data += s


class TopLevelTests(TestCase):

    """Tests of top-level functions"""

    def setUp(self):
        """ """
        self.in_seqs_variable_len_bc1 = in_seqs_variable_len_bc1
        self.bc_map_variable_len_bc1 = bc_map_variable_len_bc1
        self.primer_seq_lens_variable_len_bc1 = primer_seq_lens_variable_len_bc1
        self.all_primers_variable_len_bc1 = all_primers_variable_len_bc1
        self.expected_fasta_variable_len_bc1 = expected_fasta_variable_len_bc1
        self.in_seqs_variable_len_bc2 = in_seqs_variable_len_bc1
        self.bc_map_variable_len_bc2 = bc_map_variable_len_bc2
        self.primer_seq_lens_variable_len_bc2 = primer_seq_lens_variable_len_bc2
        self.all_primers_variable_len_bc2 = all_primers_variable_len_bc2
        self.expected_fasta_variable_len_bc2 = expected_fasta_variable_len_bc1
        self.in_seqs_fixed_len_bc1 = in_seqs_fixed_len_bc1
        self.bc_map_fixed_len_bc1 = bc_map_fixed_len_bc1
        self.primer_seq_lens_fixed_len_bc1 = primer_seq_lens_fixed_len_bc1
        self.all_primers_fixed_len_bc1 = all_primers_fixed_len_bc1
        self.expected_fasta_fixed_len_bc1 = expected_fasta_fixed_len_bc1
        self.in_seqs_fixed_len_bc2 = in_seqs_fixed_len_bc1
        self.bc_map_fixed_len_bc2 = bc_map_fixed_len_bc1
        self.primer_seq_lens_fixed_len_bc2 = primer_seq_lens_fixed_len_bc1
        self.all_primers_fixed_len_bc2 = all_primers_fixed_len_bc1
        self.expected_fasta_fixed_len_bc2 = expected_fasta_fixed_len_bc2
        self.expected_fasta_fixed_len_bc1_no_primers = \
            expected_fasta_fixed_len_bc1_no_primers
        self.reverse_primers_fixed_len_bc1 =\
            reverse_primers_fixed_len_bc1
        self.in_seqs_reverse_primers =\
            in_seqs_reverse_primers
        self.expected_in_seqs_reverse_primers =\
            expected_in_seqs_reverse_primers
        self.in_seqs_reverse_primers_mismatch =\
            in_seqs_reverse_primers_mismatch
        self.expected_in_seqs_reverse_primers_mismatch =\
            expected_in_seqs_reverse_primers_mismatch
        self.expected_in_seqs_reverse_primers_full_remove =\
            expected_in_seqs_reverse_primers_full_remove
        self.expected_in_seqs_reverse_primers_mismatch_allowed =\
            expected_in_seqs_reverse_primers_mismatch_allowed
        self.expected_fasta_fixed_len_bc1_sliding_window =\
            expected_fasta_fixed_len_bc1_sliding_window
        self.in_seqs_fixed_len_extra_bc = in_seqs_fixed_len_extra_bc
        self.expected_fasta_extra_bc = expected_fasta_extra_bc
        self.expected_fasta_mad = expected_fasta_mad
        self.sample_seqs_fna_file = sample_seqs_fna_file
        self.expected_fasta_fixed_added_demultiplex =\
            expected_fasta_fixed_added_demultiplex
        self.in_seqs_added_demultiplex = in_seqs_added_demultiplex
        self.bc_map_added_demultiplex = bc_map_added_demultiplex
        self.bc_map_added_demultiplex_group = bc_map_added_demultiplex_group
        self.expected_fasta_added_demultiplex_group =\
            expected_fasta_added_demultiplex_group
        self.in_seqs_fixed_len_bc1_qual_scores =\
            in_seqs_fixed_len_bc1_qual_scores
        self.sample_mapping = sample_mapping
        self.sample_mapping_var_length = sample_mapping_var_length
        self.sample_mapping_bad_char_sampleid = sample_mapping_bad_char_sampleid
        self.sample_mapping_bad_char_datafield =\
            sample_mapping_bad_char_datafield
        self.in_seqs_ambi_chars = in_seqs_ambi_chars

        fd, self.sample_fasta_file = mkstemp(prefix="sample_seqs_",
                                            suffix=".fasta")
        close(fd)
        seq_file = open(self.sample_fasta_file, 'w')
        seq_file.write("\n".join(self.in_seqs_fixed_len_extra_bc))
        seq_file.close()

        fd, self.sample_qual_file = mkstemp(prefix="sample_qual_",
                                           suffix=".qual")
        close(fd)
        qual_file = open(self.sample_qual_file, "w")
        qual_file.write("\n".join(self.in_seqs_fixed_len_bc1_qual_scores))
        qual_file.close()

        fd, self.sample_mapping_file = mkstemp(prefix="sample_mapping_",
                                              suffix=".txt")
        close(fd)
        map_file = open(self.sample_mapping_file, "w")
        map_file.write(self.sample_mapping)
        map_file.close()

        fd, self.sample_mapping_file_var_length =\
            mkstemp(
                prefix="sample_mapping_var_len",
                suffix=".txt")
        close(fd)
        map_file = open(self.sample_mapping_file_var_length, "w")
        map_file.write(self.sample_mapping_var_length)
        map_file.close()

        fd, self.sample_mapping_bad_char_sampleid_f =\
            mkstemp(prefix="sample_mapping_bad_char_sampleid_",
                    suffix=".txt")
        close(fd)
        map_file = open(self.sample_mapping_bad_char_sampleid_f, "w")
        map_file.write(self.sample_mapping_bad_char_sampleid)
        map_file.close()

        fd, self.sample_mapping_bad_char_datafield_f =\
            mkstemp(prefix="sample_mapping_bad_char_sampleid_",
                    suffix=".txt")
        close(fd)
        map_file = open(self.sample_mapping_bad_char_datafield_f, "w")
        map_file.write(self.sample_mapping_bad_char_datafield)
        map_file.close()

        fd, self.sample_fasta_ambi_chars_f =\
            mkstemp(prefix="sample_fasta_ambi_chars_",
                    suffix=".fna")
        close(fd)
        fna_file = open(self.sample_fasta_ambi_chars_f, "w")
        fna_file.write("\n".join(self.in_seqs_ambi_chars))
        fna_file.close()

        self.output_dir = mkdtemp(prefix="split_libraries_",
                                     suffix="/")

        self._files_to_remove = \
            [self.sample_fasta_file, self.sample_qual_file,
             self.sample_mapping_file, self.sample_mapping_file_var_length,
             self.sample_mapping_bad_char_sampleid_f,
             self.sample_mapping_bad_char_datafield_f,
             self.sample_fasta_ambi_chars_f]

    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)

    def test_check_window_qual_scores(self):
        """check_window_qual_scores returns False, index if window below qual
        threshold."""
        scores1 = [8, 8, 8, 8, 8, 8, 8, 2, 2, 2, 2, 2]
        self.assertEqual(check_window_qual_scores(scores1, 5, 5), (False, 5))
        self.assertEqual(check_window_qual_scores(scores1, 10, 5), (True, 2))
        # windowsize larger than qual score list works
        self.assertEqual(check_window_qual_scores(scores1, 100, 5), (True, 0))
        self.assertEqual(check_window_qual_scores([], 5, 1), True)
        # check each base  in its own window
        self.assertEqual(check_window_qual_scores(scores1, 1, 2), (True, 11))
        self.assertEqual(check_window_qual_scores(scores1, 1, 5), (False, 7))

    def test_expand_degeneracies(self):
        """expand_degeneracies should make possible strings"""
        # No expansion.
        self.assertEqual(expand_degeneracies(['ACG']), ['ACG'])

        # Expansion, single sequence.
        self.assertEqual(sorted(expand_degeneracies(['RGY'])),
                         ['AGC', 'AGT', 'GGC', 'GGT'])

        # Multiple sequences.
        self.assertEqual(sorted(expand_degeneracies(['ACGW', 'KAT'])),
                         ['ACGA', 'ACGT', 'GAT', 'TAT'])

    def test_get_infile(self):
        """get_infile should return filehandle"""
        pass  # not practically testable, but obvious file I/O

    def test_count_mismatches(self):
        """count_mismatches should count mismatches correctly"""
        self.assertEqual(count_mismatches('GG', 'GG', 10), 0)
        self.assertEqual(count_mismatches('GGG', 'AAA', 10), 3)
        self.assertEqual(count_mismatches('GGG', 'AAA', 1), 2)

    def test_ok_mm_primer(self):
        """ok_mm_primer should test that primer is within max mismatches"""
        primers = ['AAAA', 'GGGG']
        self.assertEqual(ok_mm_primer('AAAA', primers, 0), True)
        self.assertEqual(ok_mm_primer('AAAA', primers, 3), True)
        self.assertEqual(ok_mm_primer('CCCC', primers, 0), False)
        self.assertEqual(ok_mm_primer('CCCA', primers, 3), True)
        self.assertEqual(ok_mm_primer('CCCA', primers, 2), False)
        self.assertEqual(ok_mm_primer('CCGG', primers, 2), True)
        self.assertEqual(ok_mm_primer('CCGA', primers, 2), False)

    def test_check_map(self):
        """check_map should return valid barcodes as expected"""
        s = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tX\tDescription
#fake data
x\tAA\tAC\t3\tsample_x
y\t"AC"\tAC\t4\t"sample_y"
z\tGG\tGC\t5\tsample_z"""
        f = StringIO(s)
        f.name = 'test.xls'
        headers, id_map, barcode_to_sample_id, warnings, errors, \
            primer_seqs_lens, all_primers = check_map(f,
                                                      disable_primer_check=False)

        self.assertEqual(
            barcode_to_sample_id,
            {'AA': 'x',
             'AC': 'y',
             'GG': 'z'})

        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

    def test_check_map_primer_pool(self):
        """check_map should handle primer pools as expected"""
        s = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tX\tDescription
#fake data
x\tAA\tAC\t3\tsample_x
y\t"AC"\tAT,DC\t4\t"sample_y"
z\tGG\tGC\t5\tsample_z"""
        f = StringIO(s)
        f.name = 'test.xls'
        headers, id_map, barcode_to_sample_id, warnings, errors, \
            primer_seqs_lens, all_primers = check_map(f,
                                                      disable_primer_check=False)

        self.assertEqual(
            barcode_to_sample_id,
            {'AA': 'x',
             'AC': 'y',
             'GG': 'z'})
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

        # Returns all possible primers with lengths associated.
        expected_all_primers = {'AC': 2, 'GC': 2, 'AT': 2, 'TC': 2}
        self.assertEqual(all_primers, expected_all_primers)

        # Returns all primers associated with each barcode.
        expected_primer_seqs_lens = {'AA': {'AC': 2}, 'GG': {'GC': 2},
                                     'AC': {'AC': 2, 'GC': 2, 'AT': 2, 'TC': 2}}

        self.assertEqual(primer_seqs_lens, expected_primer_seqs_lens)

    def test_fasta_ids(self):
        """fasta_ids should return list of ids in fasta files, no dups"""
        first = StringIO('>x\nACT\n>y\nAAA')
        first_copy = StringIO('>x\nACT\n>y\nAAA')
        second = StringIO('>a\nGGG\n>b\nCCC')
        self.assertEqual(fasta_ids([first, second]), set(['x', 'y', 'a', 'b']))
        first.seek(0)  # need to reset so we can read it again
        self.assertRaises(ValueError, fasta_ids, [first, first_copy])

    def test_count_ambig(self):
        """count_ambig should count ambiguous bases in seq"""
        s = 'ACC'
        s2 = 'RNY'
        s3 = 'NA'
        self.assertEqual(count_ambig(s), 0)
        self.assertEqual(count_ambig(s2), 3)
        self.assertEqual(count_ambig(s3), 1)
        self.assertEqual(count_ambig(''), 0)

    def test_split_seq(self):
        """split_seq should split seq into pieces"""
        seq = 'AAAACCCCCGTGTGTGT'
        barcode, primer, remainder = split_seq(seq, 4, 5)
        self.assertEqual(barcode, 'AAAA')
        self.assertEqual(primer, 'CCCCC')
        self.assertEqual(remainder, 'GTGTGTGT')

    def test_primer_exceeds_mismatches(self):
        """primer_exceeds_mismatches returns True if too many mismatches"""
        primers = ['AAAA', 'TTTT']
        exact = 'AAAA'
        mismatch_ok = 'AAAT'
        mismatch_bad = 'GGGG'
        self.assertEqual(primer_exceeds_mismatches(exact, primers, 0), False)
        self.assertEqual(primer_exceeds_mismatches(mismatch_ok, primers, 1),
                         False)
        self.assertEqual(primer_exceeds_mismatches(mismatch_bad, primers, 2),
                         True)

    # Tests for local alignment functions
    def test_local_align_primer_seq_fwd_rev_match(self):
        "local_align function can handle fwd/rev primers with no mismatches"
        # forward primer
        primer = 'TAGC'
        seq = 'TAGC'
        # mismatch_count, hit_start
        expected = (0, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        primer = 'TAGC'
        seq = 'TAGCCCCC'
        # mismatch_count, hit_start
        expected = (0, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        primer = 'TAGC'
        seq = 'CCCTAGCCCCC'
        # mismatch_count, hit_start
        expected = (0, 3)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        # different length primer
        primer = 'GTTTAGC'
        seq = 'GTTTAGC'
        # mismatch_count, hit_start
        expected = (0, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        primer = 'GCTC'
        seq = 'TAGCCCCC'
        # mismatch_count, hit_start
        expected = (1, 2)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        primer = 'GCTA'
        seq = 'CCCTAGCCCCC'
        # mismatch_count, hit_start
        expected = (1, 1)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

    def test_local_align_primer_seq_fwd_rev_match_ambig(self):
        "local_align function can handle fwd/rev primers with ambigs"
        primer = 'TASC'
        seq = 'TAGC'
        # primer_hit, target, mismatch_count, hit_start
        expected = (0, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

    def test_local_align_primer_seq_mm(self):
        "local_align function can handle fwd/rev primers with mismatches"
        # forward primer
        primer = 'AAAAACTTTTT'
        seq = 'AAAAAGTTTTT'
        # mismatch_count, hit_start
        expected = (1, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        # forward primer
        primer = 'AAAACCTTTTT'
        seq = 'AAAAAGTTTTT'
        # mismatch_count, hit_start
        expected = (2, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

    def test_local_align_primer_seq_indels_middle(self):
        "local_align function can handle fwd/rev primers with indels in middle of seq"

        # Insertion in target sequence
        primer = 'CGAATCGCTATCG'
        seq = 'CGAATCTGCTATCG'
        # mismatch count, hit_start
        expected = (1, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        # Deletion in target sequence
        primer = 'CGAATCGCTATCG'
        seq = 'CGAATGCTATCG'
        # mismatch_count, hit_start
        expected = (1, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

    def test_local_align_primer_seq_multiple_mismatch_indel(self):
        "local_align function can handle fwd/rev primers with indels and mismatches"
        # multiple insertions
        primer = 'ATCGGGCGATCATT'
        seq = 'ATCGGGTTCGATCATT'
        # mismatch_count, hit_start
        expected = (2, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        # two deletions
        primer = 'ACGGTACAGTGG'
        seq = 'ACGGCAGTGG'
        # mismatch_count, hit_start
        expected = (2, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        # deletion and mismatch
        primer = 'CATCGTCGATCA'
        seq = 'CCTCGTGATCA'
        # mismatch_count, hit_start
        expected = (2, 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

    def test_seq_exceeds_homopolymers(self):
        """seq_exceeds_homopolymers returns True if too many homopolymers"""
        self.assertEqual(seq_exceeds_homopolymers('AAACGA', 3), False)
        self.assertEqual(seq_exceeds_homopolymers('AAACGA', 2), True)
        self.assertEqual(seq_exceeds_homopolymers('AAACGA', 1), True)
        self.assertEqual(seq_exceeds_homopolymers('AAACGATTTT', 3), True)

    def test_check_barcode(self):
        """check_barcode should return False if barcode ok, True otherwise"""
        self.assertEqual(check_barcode('AA', None, ['AA']), (False, 'AA',
                                                             False))
        self.assertEqual(check_barcode('GCATCGTCCACA', 'golay_12',
                                       ['GCATCGTCAACA']), (2, 'GCATCGTCAACA', True))
        # num errors for golay code is currently in bits
        self.assertEqual(check_barcode('GGTT', 4, ['TTTT']), (2, 'TTTT', True))

    def test_make_histograms(self):
        """make_histograms should make correct histograms"""
        raw_lengths = [90, 100, 110, 110, 130, 135]
        pre_lengths = [100, 110, 105, 130, 135]
        post_lengths = [130, 135]
        raw_hist, pre_hist, post_hist, bin_edges = \
            make_histograms(raw_lengths, pre_lengths, post_lengths)
        assert_almost_equal(pre_hist, array([0, 2, 1, 0, 2]))
        assert_almost_equal(post_hist, array([0, 0, 0, 0, 2]))
        assert_almost_equal(bin_edges, array([90, 100, 110, 120, 130, 140]))

    def test_check_seqs_sliding_window(self):
        """check_seqs handles sliding window truncations/removal """

        in_seqs = self.in_seqs_fixed_len_bc1
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_fasta_fixed_len_bc1_sliding_window

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings=parse_qual_score(
                self.in_seqs_fixed_len_bc1_qual_scores),
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False,
            qual_score_window=5,
            discard_bad_windows=False,
            min_qual_score=25,
            min_seq_len=200)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

    def test_check_seqs_variable_len_bc(self):
        """check_seqs handles variable length barcodes """

        # Simple test with variable length primers
        in_seqs = self.in_seqs_variable_len_bc1
        bc_map = self.bc_map_variable_len_bc1
        primer_seq_lens = self.primer_seq_lens_variable_len_bc1
        all_primers = self.all_primers_variable_len_bc1
        expected = self.expected_fasta_variable_len_bc1

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=None,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="variable_length",
            max_bc_errors=0,
            retain_unassigned_reads=False,
            attempt_bc_correction=False,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

        # Second test, includes truncated form of one of the barcodes-the
        # longest barcode should be found first
        in_seqs = self.in_seqs_variable_len_bc2
        bc_map = self.bc_map_variable_len_bc2
        primer_seq_lens = self.primer_seq_lens_variable_len_bc2
        all_primers = self.all_primers_variable_len_bc2
        expected = self.expected_fasta_variable_len_bc2

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=None,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="variable_length",
            max_bc_errors=0,
            retain_unassigned_reads=False,
            attempt_bc_correction=False,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

    def test_check_seqs_fixed_len_bc(self):
        """check_seqs handles fixed length barcodes """

        # Third test, fixed length barcodes, fixed length primers + one
        # degenerate test.  Should correct one of the passed barcodes
        in_seqs = self.in_seqs_fixed_len_bc1
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_fasta_fixed_len_bc1

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

        # Fourth test-set max_bc_errors to 0, and allow some primer mismatches
        in_seqs = self.in_seqs_fixed_len_bc2
        bc_map = self.bc_map_fixed_len_bc2
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc2
        all_primers = self.all_primers_fixed_len_bc2
        expected = self.expected_fasta_fixed_len_bc2

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=0.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=1,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

    def test_check_seqs_no_primers(self):
        """check_seqs handles disabled primers """

        # Fifth test, no primers, fixed length barcodes
        # Should correct one of the passed barcodes
        in_seqs = self.in_seqs_fixed_len_bc1
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = {}
        all_primers = {}
        expected = self.expected_fasta_fixed_len_bc1_no_primers

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=True,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

    def test_check_seqs_reverse_primers(self):
        """check_seqs handles truncating reverse primers """

        # Initial test, should truncate all seqs
        in_seqs = self.in_seqs_reverse_primers
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers
        rev_primers_test = self.reverse_primers_fixed_len_bc1

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='truncate_only',
            rev_primers=rev_primers_test,
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

        # Second test with a mismatch in seq a, should not find reverse primer
        # and will write out entire sequence.

        in_seqs = self.in_seqs_reverse_primers_mismatch
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers_mismatch
        rev_primers_test = self.reverse_primers_fixed_len_bc1

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='truncate_only',
            rev_primers=rev_primers_test,
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

        # With reverse_primer_mismatches allowed set to 1,
        # should restore truncation.
        in_seqs = self.in_seqs_reverse_primers_mismatch
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers_mismatch_allowed
        rev_primers_test = self.reverse_primers_fixed_len_bc1

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='truncate_only',
            rev_primers=rev_primers_test,
            qual_out=False,
            reverse_primer_mismatches=1)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

        # Testing truncate_remove, which should not write sequences where
        # the reverse primer is not found
        in_seqs = self.in_seqs_reverse_primers
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers_full_remove
        rev_primers_test = self.reverse_primers_fixed_len_bc1

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='truncate_remove',
            rev_primers=rev_primers_test,
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

        # Testing truncate_remove, with reverse_primer_mismatches set to 1
        # should allow all 4 sequences to be written, truncated
        in_seqs = self.in_seqs_reverse_primers_mismatch
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers_mismatch_allowed
        rev_primers_test = self.reverse_primers_fixed_len_bc1

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=1,
            disable_primer_check=False,
            reverse_primers='truncate_remove',
            rev_primers=rev_primers_test,
            qual_out=False,
            reverse_primer_mismatches=1)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

    def test_check_seqs_qual_out(self):
        """ check_seqs handles optional quality output file """

        in_seqs = self.in_seqs_fixed_len_bc1
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = expected_qual_fixed_len_bc1

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        qual_out_f = FakeOutFile()

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings=parse_qual_score(
                self.in_seqs_fixed_len_bc1_qual_scores),
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=qual_out_f)

        self.assertEqual(qual_out_f.data, expected)

    def test_check_seqs_retain_unassigned_reads(self):
        """ check_seqs handles retaining Unassigned reads """

        in_seqs = self.in_seqs_fixed_len_extra_bc
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_fasta_extra_bc

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=True,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

    def test_check_seqs_median_abs_dev(self):
        """ check_seqs handles median absolute deviation calculations """

        in_seqs = self.in_seqs_fixed_len_bc1
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_fasta_mad

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False,
            median_length_filtering=2.0)

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

    def test_check_seqs_added_demultiplex(self):
        """check_seqs handles added demultiplex field"""

        # Test added demultiplex for the run_prefix
        in_seqs = self.in_seqs_added_demultiplex
        bc_map = self.bc_map_added_demultiplex
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_fasta_fixed_added_demultiplex

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False,
            added_demultiplex_field='run_prefix')

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

        # Demultiplex by the 'group' in the fasta label
        in_seqs = self.in_seqs_added_demultiplex
        bc_map = self.bc_map_added_demultiplex_group
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_fasta_added_demultiplex_group

        fd, out_fp = mkstemp(prefix="sample_seqs_", suffix=".fna.tmp")
        close(fd)
        out_f = open(out_fp, "w")
        self._files_to_remove.append(out_f.name.replace('.tmp', ''))

        actual = check_seqs(
            fasta_out=out_f,
            fasta_files=[in_seqs],
            starting_ix=0,
            valid_map=bc_map,
            qual_mappings={},
            filters=[],
            barcode_len=12,
            keep_primer=False,
            keep_barcode=False,
            barcode_type="golay_12",
            max_bc_errors=1.5,
            retain_unassigned_reads=False,
            attempt_bc_correction=True,
            primer_seqs_lens=primer_seq_lens,
            all_primers=all_primers,
            max_primer_mm=0,
            disable_primer_check=False,
            reverse_primers='disable',
            rev_primers={},
            qual_out=False,
            added_demultiplex_field='group')

        out_f = open(out_f.name.replace('.tmp', ''), "U")
        actual_results = '\n'.join([line.strip() for line in out_f])

        self.assertEqual(actual_results, expected)

    def test_preprocess_variable_length_barcodes(self):
        """ Tests overall module with variable length barcodes """

         # Should discard all reads due to sequence length being too short

        fasta_files = [self.sample_fasta_file]
        qual_files = [self.sample_qual_file]
        mapping_file = self.sample_mapping_file_var_length
        barcode_type = "variable_length"
        min_seq_len = 200
        max_seq_len = 1000
        min_qual_score = 25
        starting_ix = 1
        keep_primer = False
        max_ambig = 0
        max_primer_mm = 1
        trim_seq_len = True
        dir_prefix = self.output_dir
        max_bc_errors = 2
        max_homopolymer = 4
        retain_unassigned_reads = False
        keep_barcode = False
        attempt_bc_correction = True
        qual_score_window = 0
        disable_primer_check = False
        reverse_primers = 'disable'
        record_qual_scores = False
        discard_bad_windows = False
        median_length_filtering = None
        added_demultiplex_field = None

        preprocess(fasta_files,
                   qual_files,
                   mapping_file,
                   barcode_type,
                   min_seq_len,
                   max_seq_len,
                   min_qual_score,
                   starting_ix,
                   keep_primer,
                   max_ambig,
                   max_primer_mm,
                   trim_seq_len,
                   dir_prefix,
                   max_bc_errors,
                   max_homopolymer,
                   retain_unassigned_reads,
                   keep_barcode,
                   attempt_bc_correction,
                   qual_score_window,
                   disable_primer_check,
                   reverse_primers,
                   record_qual_scores,
                   discard_bad_windows,
                   median_length_filtering,
                   added_demultiplex_field)

        output_seqs = open(dir_prefix + "seqs.fna", "U")
        output_log = open(dir_prefix + "split_library_log.txt", "U")
        output_histograms = open(dir_prefix + "histograms.txt", "U")

        actual_seqs = [line for line in output_seqs]
        actual_log = [line for line in output_log]
        actual_histograms = [line for line in output_histograms]

        expected_seqs = []
        expected_log = [
            'Number raw input seqs\t6\n',
            '\n',
            'Length outside bounds of 200 and 1000\t6\n',
            'Num ambiguous bases exceeds limit of 0\t0\n',
            'Missing Qual Score\t0\n',
            'Mean qual score below minimum of 25\t0\n',
            'Max homopolymer run exceeds limit of 4\t0\n',
            'Num mismatches in primer exceeds limit of 1: 0\n',
            '\n',
            'Sequence length details for all sequences passing quality filters:\n',
            'No sequences passed quality filters for writing.\n',
            '\n',
            'Barcodes corrected/not\t0/0\n',
            'Uncorrected barcodes will not be written to the output fasta file.\n',
            'Corrected barcodes will be written with the appropriate barcode category.\n',
            'Corrected but unassigned sequences will not be written unless --retain_unassigned_reads is enabled.\n',
            '\n',
            'Total valid barcodes that are not in mapping file\t0\n',
            'Sequences associated with valid barcodes that are not in the mapping file will not be written.\n',
            '\n',
            'Barcodes in mapping file\n',
            'Sample\tSequence Count\tBarcode\n',
            's2\t0\tAGAGTCCTGAGC\n',
            's1\t0\tACACATGTCTA\n',
            's3\t0\tAACTGTGCGTACG\n',
            '\n',
            'Total number seqs written\t0']
        expected_histograms = [
            '# bins raw sequence lengths, length of sequences that pass quality filters before processing, and lengths of sequences that pass quality filters post processing.\n',
            'Length\tRaw\tBefore\tAfter\n',
            '20\t2\t0\t0\n',
            '30\t4\t0\t0']

        self.assertEqual(actual_seqs, expected_seqs)
        self.assertEqual(actual_log, expected_log)
        self.assertEqual(actual_histograms, expected_histograms)

    def test_preprocess(self):
        """ Overall module functionality test """

        # Should discard all reads due to sequence length being too short

        fasta_files = [self.sample_fasta_file]
        qual_files = [self.sample_qual_file]
        mapping_file = self.sample_mapping_file
        barcode_type = "golay_12"
        min_seq_len = 200
        max_seq_len = 1000
        min_qual_score = 25
        starting_ix = 1
        keep_primer = False
        max_ambig = 0
        max_primer_mm = 1
        trim_seq_len = True
        dir_prefix = self.output_dir
        max_bc_errors = 2
        max_homopolymer = 4
        retain_unassigned_reads = False
        keep_barcode = False
        attempt_bc_correction = True
        qual_score_window = 0
        disable_primer_check = False
        reverse_primers = 'disable'
        record_qual_scores = False
        discard_bad_windows = False
        median_length_filtering = None
        added_demultiplex_field = None

        preprocess(fasta_files,
                   qual_files,
                   mapping_file,
                   barcode_type,
                   min_seq_len,
                   max_seq_len,
                   min_qual_score,
                   starting_ix,
                   keep_primer,
                   max_ambig,
                   max_primer_mm,
                   trim_seq_len,
                   dir_prefix,
                   max_bc_errors,
                   max_homopolymer,
                   retain_unassigned_reads,
                   keep_barcode,
                   attempt_bc_correction,
                   qual_score_window,
                   disable_primer_check,
                   reverse_primers,
                   record_qual_scores,
                   discard_bad_windows,
                   median_length_filtering,
                   added_demultiplex_field)

        output_seqs = open(dir_prefix + "seqs.fna", "U")
        output_log = open(dir_prefix + "split_library_log.txt", "U")
        output_histograms = open(dir_prefix + "histograms.txt", "U")

        actual_seqs = [line for line in output_seqs]
        actual_log = [line for line in output_log]
        actual_histograms = [line for line in output_histograms]

        expected_seqs = []
        expected_log = [
            'Number raw input seqs\t6\n',
            '\n',
            'Length outside bounds of 200 and 1000\t6\n',
            'Num ambiguous bases exceeds limit of 0\t0\n',
            'Missing Qual Score\t0\n',
            'Mean qual score below minimum of 25\t0\n',
            'Max homopolymer run exceeds limit of 4\t0\n',
            'Num mismatches in primer exceeds limit of 1: 0\n',
            '\n',
            'Sequence length details for all sequences passing quality filters:\n',
            'No sequences passed quality filters for writing.\n',
            '\n',
            'Barcodes corrected/not\t0/0\n',
            'Uncorrected barcodes will not be written to the output fasta file.\n',
            'Corrected barcodes will be written with the appropriate barcode category.\n',
            'Corrected but unassigned sequences will not be written unless --retain_unassigned_reads is enabled.\n',
            '\n',
            'Total valid barcodes that are not in mapping file\t0\n',
            'Sequences associated with valid barcodes that are not in the mapping file will not be written.\n',
            '\n',
            'Barcodes in mapping file\n',
            'Sample\tSequence Count\tBarcode\n',
            's2\t0\tAGAGTCCTGAGC\n',
            's1\t0\tACACATGTCTAC\n',
            's3\t0\tAACTGTGCGTAC\n',
            '\n',
            'Total number seqs written\t0']
        expected_histograms = [
            '# bins raw sequence lengths, length of sequences that pass quality filters before processing, and lengths of sequences that pass quality filters post processing.\n',
            'Length\tRaw\tBefore\tAfter\n',
            '20\t2\t0\t0\n',
            '30\t4\t0\t0']

        self.assertEqual(actual_seqs, expected_seqs)
        self.assertEqual(actual_log, expected_log)
        self.assertEqual(actual_histograms, expected_histograms)

        # With minimal length at 5, should retain 4 sequences

        fasta_files = [self.sample_fasta_file]
        qual_files = [self.sample_qual_file]
        mapping_file = self.sample_mapping_file
        barcode_type = "golay_12"
        min_seq_len = 5
        max_seq_len = 1000
        min_qual_score = 25
        starting_ix = 1
        keep_primer = False
        max_ambig = 0
        max_primer_mm = 0
        trim_seq_len = False
        dir_prefix = self.output_dir
        max_bc_errors = 2
        max_homopolymer = 4
        retain_unassigned_reads = False
        keep_barcode = False
        attempt_bc_correction = True
        qual_score_window = 0
        disable_primer_check = False
        reverse_primers = 'disable'
        record_qual_scores = False
        discard_bad_windows = False
        median_length_filtering = None
        added_demultiplex_field = None

        preprocess(fasta_files,
                   qual_files,
                   mapping_file,
                   barcode_type,
                   min_seq_len,
                   max_seq_len,
                   min_qual_score,
                   starting_ix,
                   keep_primer,
                   max_ambig,
                   max_primer_mm,
                   trim_seq_len,
                   dir_prefix,
                   max_bc_errors,
                   max_homopolymer,
                   retain_unassigned_reads,
                   keep_barcode,
                   attempt_bc_correction,
                   qual_score_window,
                   disable_primer_check,
                   reverse_primers,
                   record_qual_scores,
                   discard_bad_windows,
                   median_length_filtering,
                   added_demultiplex_field)

        output_seqs = open(dir_prefix + "seqs.fna", "U")
        output_log = open(dir_prefix + "split_library_log.txt", "U")
        output_histograms = open(dir_prefix + "histograms.txt", "U")

        actual_seqs = [line for line in output_seqs]
        actual_log = [line for line in output_log]
        actual_histograms = [line for line in output_histograms]

        expected_seqs = [
            '>s1_1 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\n',
            'CCCTTATATATATAT\n',
            '>s2_2 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\n',
            'CCCTTTCCA\n',
            '>s3_3 c orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0\n',
            'AACCGGCCGGTT\n',
            '>s1_4 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\n',
            'CCCTTACTATATAT\n']
        expected_log = [
            'Number raw input seqs\t6\n',
            '\n',
            'Length outside bounds of 5 and 1000\t0\n',
            'Num ambiguous bases exceeds limit of 0\t0\n',
            'Missing Qual Score\t0\n',
            'Mean qual score below minimum of 25\t0\n',
            'Max homopolymer run exceeds limit of 4\t0\n',
            'Num mismatches in primer exceeds limit of 0: 2\n',
            '\n',
            'Sequence length details for all sequences passing quality filters:\n',
            'Raw len min/max/avg\t29.0/35.0/32.5\n',
            'Wrote len min/max/avg\t9.0/15.0/12.5\n',
            '\n',
            'Barcodes corrected/not\t1/0\n',
            'Uncorrected barcodes will not be written to the output fasta file.\n',
            'Corrected barcodes will be written with the appropriate barcode category.\n',
            'Corrected but unassigned sequences will not be written unless --retain_unassigned_reads is enabled.\n',
            '\n',
            'Total valid barcodes that are not in mapping file\t0\n',
            'Sequences associated with valid barcodes that are not in the mapping file will not be written.\n',
            '\n',
            'Barcodes in mapping file\n',
            'Num Samples\t3\n',
            'Sample ct min/max/mean: 1 / 2 / 1.33\n',
            'Sample\tSequence Count\tBarcode\n',
            's1\t2\tACACATGTCTAC\n',
            's2\t1\tAGAGTCCTGAGC\n',
            's3\t1\tAACTGTGCGTAC\n',
            '\n',
            'Total number seqs written\t4']
        expected_histograms = [
            '# bins raw sequence lengths, length of sequences that pass quality filters before processing, and lengths of sequences that pass quality filters post processing.\n',
            'Length\tRaw\tBefore\tAfter\n',
            '0\t0\t0\t1\n',
            '10\t0\t0\t3\n',
            '20\t2\t1\t0\n',
            '30\t4\t3\t0']

        self.assertEqual(actual_seqs, expected_seqs)
        self.assertEqual(actual_log, expected_log)
        self.assertEqual(actual_histograms, expected_histograms)

        # Added sliding window should discard read "b"

        fasta_files = [self.sample_fasta_file]
        qual_files = [self.sample_qual_file]
        mapping_file = self.sample_mapping_file
        barcode_type = "golay_12"
        min_seq_len = 5
        max_seq_len = 1000
        min_qual_score = 22
        starting_ix = 1
        keep_primer = False
        max_ambig = 0
        max_primer_mm = 0
        trim_seq_len = False
        dir_prefix = self.output_dir
        max_bc_errors = 2
        max_homopolymer = 4
        retain_unassigned_reads = False
        keep_barcode = False
        attempt_bc_correction = True
        qual_score_window = 3
        disable_primer_check = False
        reverse_primers = 'disable'
        record_qual_scores = False
        discard_bad_windows = True
        median_length_filtering = None
        added_demultiplex_field = None
        reverse_primer_mismatches = 0

        preprocess(fasta_files,
                   qual_files,
                   mapping_file,
                   barcode_type,
                   min_seq_len,
                   max_seq_len,
                   min_qual_score,
                   starting_ix,
                   keep_primer,
                   max_ambig,
                   max_primer_mm,
                   trim_seq_len,
                   dir_prefix,
                   max_bc_errors,
                   max_homopolymer,
                   retain_unassigned_reads,
                   keep_barcode,
                   attempt_bc_correction,
                   qual_score_window,
                   disable_primer_check,
                   reverse_primers,
                   reverse_primer_mismatches,
                   record_qual_scores,
                   discard_bad_windows,
                   median_length_filtering,
                   added_demultiplex_field)

        output_seqs = open(dir_prefix + "seqs.fna", "U")
        output_log = open(dir_prefix + "split_library_log.txt", "U")
        output_histograms = open(dir_prefix + "histograms.txt", "U")

        actual_seqs = [line for line in output_seqs]
        actual_log = [line for line in output_log]
        actual_histograms = [line for line in output_histograms]

        expected_seqs = [
            '>s1_1 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\n',
            'CCCTTATATATATAT\n',
            '>s3_2 c orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0\n',
            'AACCGGCCGGTT\n',
            '>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\n',
            'CCCTTACTATATAT\n']
        expected_log = [
            'Number raw input seqs\t6\n',
            '\n',
            'Length outside bounds of 5 and 1000\t0\n',
            'Num ambiguous bases exceeds limit of 0\t0\n',
            'Missing Qual Score\t0\n',
            'Mean qual score below minimum of 22\t0\n',
            'Max homopolymer run exceeds limit of 4\t0\n',
            'Num mismatches in primer exceeds limit of 0: 2\n',
            '\n',
            'Size of quality score window, in base pairs: 3\n',
            'Number of sequences where a low quality score window was detected: 1\n',
            'Sequences with a low quality score were not written, -g option enabled.\n',
            '\n',
            'Sequence length details for all sequences passing quality filters:\n',
            'Raw len min/max/avg\t32.0/35.0/33.7\n',
            'Wrote len min/max/avg\t12.0/15.0/13.7\n',
            '\n',
            'Barcodes corrected/not\t1/0\n',
            'Uncorrected barcodes will not be written to the output fasta file.\n',
            'Corrected barcodes will be written with the appropriate barcode category.\n',
            'Corrected but unassigned sequences will not be written unless --retain_unassigned_reads is enabled.\n',
            '\n',
            'Total valid barcodes that are not in mapping file\t0\n',
            'Sequences associated with valid barcodes that are not in the mapping file will not be written.\n',
            '\n',
            'Barcodes in mapping file\n',
            'Num Samples\t2\n',
            'Sample ct min/max/mean: 1 / 2 / 1.50\n',
            'Sample\tSequence Count\tBarcode\n',
            's1\t2\tACACATGTCTAC\n',
            's3\t1\tAACTGTGCGTAC\n',
            's2\t0\tAGAGTCCTGAGC\n',
            '\n',
            'Total number seqs written\t3']
        expected_histograms = [
            '# bins raw sequence lengths, length of sequences that pass quality filters before processing, and lengths of sequences that pass quality filters post processing.\n',
            'Length\tRaw\tBefore\tAfter\n',
            '10\t0\t0\t3\n',
            '20\t2\t0\t0\n',
            '30\t4\t3\t0']

        self.assertEqual(actual_seqs, expected_seqs)
        self.assertEqual(actual_log, expected_log)
        self.assertEqual(actual_histograms, expected_histograms)

    def test_preprocess_ambi_trunc(self):
        """ Overall module test for 'N' character truncation """

        # Will truncate sequences at the first N character, remove one seq
        # due to being less than 30 nucleotides following truncation
        # counting the sequence length of the barcodes + primer

        fasta_files = [self.sample_fasta_ambi_chars_f]
        qual_files = [self.sample_qual_file]
        mapping_file = self.sample_mapping_file
        barcode_type = "golay_12"
        min_seq_len = 30
        max_seq_len = 1000
        min_qual_score = 22
        starting_ix = 1
        keep_primer = False
        max_ambig = 0
        max_primer_mm = 0
        trim_seq_len = False
        dir_prefix = self.output_dir
        max_bc_errors = 2
        max_homopolymer = 4
        retain_unassigned_reads = False
        keep_barcode = False
        attempt_bc_correction = True
        qual_score_window = False
        disable_primer_check = False
        reverse_primers = 'disable'
        record_qual_scores = False
        discard_bad_windows = False
        median_length_filtering = None
        added_demultiplex_field = None
        reverse_primer_mismatches = 0
        truncate_ambi_bases = True

        preprocess(fasta_files,
                   qual_files,
                   mapping_file,
                   barcode_type,
                   min_seq_len,
                   max_seq_len,
                   min_qual_score,
                   starting_ix,
                   keep_primer,
                   max_ambig,
                   max_primer_mm,
                   trim_seq_len,
                   dir_prefix,
                   max_bc_errors,
                   max_homopolymer,
                   retain_unassigned_reads,
                   keep_barcode,
                   attempt_bc_correction,
                   qual_score_window,
                   disable_primer_check,
                   reverse_primers,
                   reverse_primer_mismatches,
                   record_qual_scores,
                   discard_bad_windows,
                   median_length_filtering,
                   added_demultiplex_field,
                   truncate_ambi_bases)

        output_seqs = open(dir_prefix + "seqs.fna", "U")
        output_log = open(dir_prefix + "split_library_log.txt", "U")
        output_histograms = open(dir_prefix + "histograms.txt", "U")

        actual_seqs = [line for line in output_seqs]
        actual_log = [line for line in output_log]
        actual_histograms = [line for line in output_histograms]

        expected_seqs = [
            '>s1_1 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\n',
            'CCCTTATATATAT\n',
            '>s3_2 c orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0\n',
            'AACCGGCCGGTT\n',
            '>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\n',
            'CCCTTACTACCGA\n']
        expected_log = [
            'Number raw input seqs\t6\n',
            '\n',
            'Length outside bounds of 30 and 1000\t2\n',
            'Missing Qual Score\t0\n',
            'Mean qual score below minimum of 22\t0\n',
            'Max homopolymer run exceeds limit of 4\t0\n',
            'Num mismatches in primer exceeds limit of 0: 1\n',
            '\n',
            'Truncation at first ambiguous "N" character enabled.\n',
            'Sequences discarded after truncation due to sequence length below the minimum 30: 0\n',
            '\n',
            'Sequence length details for all sequences passing quality filters:\n',
            'Raw len min/max/avg\t32.0/39.0/35.7\n',
            'Wrote len min/max/avg\t12.0/13.0/12.7\n',
            '\n',
            'Barcodes corrected/not\t1/0\n',
            'Uncorrected barcodes will not be written to the output fasta file.\n',
            'Corrected barcodes will be written with the appropriate barcode category.\n',
            'Corrected but unassigned sequences will not be written unless --retain_unassigned_reads is enabled.\n',
            '\n',
            'Total valid barcodes that are not in mapping file\t0\n',
            'Sequences associated with valid barcodes that are not in the mapping file will not be written.\n',
            '\n',
            'Barcodes in mapping file\n',
            'Num Samples\t2\n',
            'Sample ct min/max/mean: 1 / 2 / 1.50\n',
            'Sample\tSequence Count\tBarcode\n',
            's1\t2\tACACATGTCTAC\n',
            's3\t1\tAACTGTGCGTAC\n',
            's2\t0\tAGAGTCCTGAGC\n',
            '\n',
            'Total number seqs written\t3']
        expected_histograms = [
            '# bins raw sequence lengths, length of sequences that pass quality filters before processing, and lengths of sequences that pass quality filters post processing.\n',
            'Length\tRaw\tBefore\tAfter\n',
            '10\t0\t0\t3\n',
            '20\t2\t0\t0\n',
            '30\t4\t3\t0']

        self.assertEqual(actual_seqs, expected_seqs)
        self.assertEqual(actual_log, expected_log)
        self.assertEqual(actual_histograms, expected_histograms)

    def test_preprocess_bad_chars_in_mapping(self):
        """ Overall module functionality test with invalid characters """

        # Should discard all reads due to sequence length being too short
        # But should not halt due to bad characters in a data field

        fasta_files = [self.sample_fasta_file]
        qual_files = [self.sample_qual_file]
        mapping_file = self.sample_mapping_bad_char_datafield_f
        barcode_type = "golay_12"
        min_seq_len = 200
        max_seq_len = 1000
        min_qual_score = 25
        starting_ix = 1
        keep_primer = False
        max_ambig = 0
        max_primer_mm = 1
        trim_seq_len = True
        dir_prefix = self.output_dir
        max_bc_errors = 2
        max_homopolymer = 4
        retain_unassigned_reads = False
        keep_barcode = False
        attempt_bc_correction = True
        qual_score_window = 0
        disable_primer_check = False
        reverse_primers = 'disable'
        record_qual_scores = False
        discard_bad_windows = False
        median_length_filtering = None
        added_demultiplex_field = None

        preprocess(fasta_files,
                   qual_files,
                   mapping_file,
                   barcode_type,
                   min_seq_len,
                   max_seq_len,
                   min_qual_score,
                   starting_ix,
                   keep_primer,
                   max_ambig,
                   max_primer_mm,
                   trim_seq_len,
                   dir_prefix,
                   max_bc_errors,
                   max_homopolymer,
                   retain_unassigned_reads,
                   keep_barcode,
                   attempt_bc_correction,
                   qual_score_window,
                   disable_primer_check,
                   reverse_primers,
                   record_qual_scores,
                   discard_bad_windows,
                   median_length_filtering,
                   added_demultiplex_field)

        output_seqs = open(dir_prefix + "seqs.fna", "U")
        output_log = open(dir_prefix + "split_library_log.txt", "U")
        output_histograms = open(dir_prefix + "histograms.txt", "U")

        actual_seqs = [line for line in output_seqs]
        actual_log = [line for line in output_log]
        actual_histograms = [line for line in output_histograms]

        expected_seqs = []
        expected_log = [
            'Number raw input seqs\t6\n',
            '\n',
            'Length outside bounds of 200 and 1000\t6\n',
            'Num ambiguous bases exceeds limit of 0\t0\n',
            'Missing Qual Score\t0\n',
            'Mean qual score below minimum of 25\t0\n',
            'Max homopolymer run exceeds limit of 4\t0\n',
            'Num mismatches in primer exceeds limit of 1: 0\n',
            '\n',
            'Sequence length details for all sequences passing quality filters:\n',
            'No sequences passed quality filters for writing.\n',
            '\n',
            'Barcodes corrected/not\t0/0\n',
            'Uncorrected barcodes will not be written to the output fasta file.\n',
            'Corrected barcodes will be written with the appropriate barcode category.\n',
            'Corrected but unassigned sequences will not be written unless --retain_unassigned_reads is enabled.\n',
            '\n',
            'Total valid barcodes that are not in mapping file\t0\n',
            'Sequences associated with valid barcodes that are not in the mapping file will not be written.\n',
            '\n',
            'Barcodes in mapping file\n',
            'Sample\tSequence Count\tBarcode\n',
            's2\t0\tAGAGTCCTGAGC\n',
            's1\t0\tACACATGTCTAC\n',
            's3\t0\tAACTGTGCGTAC\n',
            '\n',
            'Total number seqs written\t0']
        expected_histograms = [
            '# bins raw sequence lengths, length of sequences that pass quality filters before processing, and lengths of sequences that pass quality filters post processing.\n',
            'Length\tRaw\tBefore\tAfter\n',
            '20\t2\t0\t0\n',
            '30\t4\t0\t0']

        self.assertEqual(actual_seqs, expected_seqs)
        self.assertEqual(actual_log, expected_log)
        self.assertEqual(actual_histograms, expected_histograms)

        '''# With invalid character in a SampleID, should raise ValueError

        fasta_files = [self.sample_fasta_file]
        qual_files = [self.sample_qual_file]
        mapping_file = self.sample_mapping_bad_char_sampleid_f
        barcode_type="golay_12"
        min_seq_len=200
        max_seq_len=1000
        min_qual_score=25
        starting_ix=1
        keep_primer=False
        max_ambig=0
        max_primer_mm=1
        trim_seq_len=True
        dir_prefix=self.output_dir
        max_bc_errors=2
        max_homopolymer=4
        retain_unassigned_reads=False
        keep_barcode=False
        attempt_bc_correction=True
        qual_score_window=0
        disable_primer_check=False
        reverse_primers='disable'
        record_qual_scores=False
        discard_bad_windows=False
        median_length_filtering=None
        added_demultiplex_field=None


        self.assertRaises(ValueError, preprocess, fasta_files,
                           qual_files,
                           mapping_file,
                           barcode_type,
                           min_seq_len,
                           max_seq_len,
                           min_qual_score,
                           starting_ix,
                           keep_primer,
                           max_ambig,
                           max_primer_mm,
                           trim_seq_len,
                           dir_prefix,
                           max_bc_errors,
                           max_homopolymer,
                           retain_unassigned_reads,
                           keep_barcode,
                           attempt_bc_correction,
                           qual_score_window,
                           disable_primer_check,
                           reverse_primers,
                           record_qual_scores,
                           discard_bad_windows,
                           median_length_filtering,
                           added_demultiplex_field)'''


in_seqs_variable_len_bc1 = """>a
ACCGGTCCGGACCCTTATATATATAT
>b
AGGAGTCCGGACCCTTTCCA
>c
ATTAACCCGGAAACCGGCCGGTT
>d
ACCGGTCCGGACCCTTACTATATAT
>e
TTTTGTCCGGACCCTTACTATATAT
>d_primer_error
ACCTGGTCCGGACCCTTACTATATAT
""".split()

bc_map_variable_len_bc1 = {'ACC': 's1', 'AGGA': 's2', 'ATTA': 's3'}
primer_seq_lens_variable_len_bc1 = {'ACC': {'GGTCCGGA': 8},
                                    'AGGA': {'GTCCGGA': 7},
                                    'ATTA': {'ACCCGGA': 7}}
all_primers_variable_len_bc1 = {'GGTCCGGA': 8, 'GTCCGGA': 7, 'ACCCGGA': 7}

expected_fasta_variable_len_bc1 = """>s1_0 a orig_bc=ACC new_bc=ACC bc_diffs=0
CCCTTATATATATAT
>s2_1 b orig_bc=AGGA new_bc=AGGA bc_diffs=0
CCCTTTCCA
>s3_2 c orig_bc=ATTA new_bc=ATTA bc_diffs=0
AACCGGCCGGTT
>s1_3 d orig_bc=ACC new_bc=ACC bc_diffs=0
CCCTTACTATATAT"""

bc_map_variable_len_bc2 = {
    'ACC': 's1',
    'AGGA': 's2',
    'ATTA': 's3',
    'AGG': 's4'}
primer_seq_lens_variable_len_bc2 = {'ACC': {'GGTCCGGA': 8},
                                    'AGGA': {'GTCCGGA': 7},
                                    'ATTA': {'ACCCGGA': 7},
                                    'AGG': {'AGTCCGGA': 8}}
all_primers_variable_len_bc2 = {'GGTCCGGA': 8, 'GTCCGGA': 7, 'ACCCGGA': 7,
                                'AGTCCGGA': 8}

# Fixed barcode test data
in_seqs_fixed_len_bc1 = """>a
ACACATGTCTACGGTCCGGACCCTTATATATATAT
>b
AGAGTCCTGAGCGGTCCGGACCCTTTCCA
>c
AATCGTGACTCGGGTCTGGAAACCGGCCGGTT
>d
ACTCATGTCTACGGTCCGGACCCTTACTATATAT
>e_no_barcode_match
TTTTGTCCGGACCCTTACTATATAT
>d_primer_error
AGAGTCCTGAGCGGTCCGGTACGTTTACTGGA
""".split('\n')

# Fixed barcode test data with "N" character included
in_seqs_ambi_chars = """>a
ACACATGTCTACGGTCCGGACCCTTATATATATNAT
>b
AGAGTCCTGAGCGGTCCGGACCCTTTCCA
>c
AACTGTGCGTACGGTCTGGAAACCGGCCGGTT
>d
ACTCATGTCTACGGTCCGGACCCTTACTACCGANTATAT
>e_no_barcode_match
TTTTGTCCGGACCCTTACTATATAT
>d_primer_error
AGAGTCCTGAGCGGTCCGGTACGTTTACTGGA
""".split('\n')

sample_mapping = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDescription
# Test mapping for split_libraries.py unit tests
s1\tACACATGTCTAC\tGGTCCGGA\tControl\ts1_mouse
s2\tAGAGTCCTGAGC\tGGTCCGGA\tControl\ts2_mouse
s3\tAACTGTGCGTAC\tGGTCYGGA\tFasted\ts3_mouse
"""

# Sample mapping with invalid character in one sampleID
sample_mapping_bad_char_sampleid = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDescription
# Test mapping for split_libraries.py unit tests
s1\tACACATGTCTAC\tGGTCCGGA\tControl\ts1_mouse
s_2\tAGAGTCCTGAGC\tGGTCCGGA\tControl\ts2_mouse
s3\tAACTGTGCGTAC\tGGTCYGGA\tFasted\ts3_mouse
"""

# Sample mapping with invalid character in data field
sample_mapping_bad_char_datafield = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDescription
# Test mapping for split_libraries.py unit tests
s1\tACACATGTCTAC\tGGTCCGGA\tControl\ts1_mouse
s2\tAGAGTCCTGAGC\tGGTCCGGA\tControl\ts2_mouse
s3\tAACTGTGCGTAC\tGGTCYGGA\tFas^^ted\ts3_mouse
"""

sample_mapping_var_length = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDescription
# Test mapping for split_libraries.py unit tests
s1\tACACATGTCTA\tCGGTCCGGA\tControl\ts1_mouse
s2\tAGAGTCCTGAGC\tGGTCCGGA\tControl\ts2_mouse
s3\tAACTGTGCGTACG\tGTCYGGA\tFasted\ts3_mouse
"""


# Fixed barcode test data, with one valid BC not in mapping data
in_seqs_fixed_len_extra_bc = """>a
ACACATGTCTACGGTCCGGACCCTTATATATATAT
>b
AGAGTCCTGAGCGGTCCGGACCCTTTCCA
>c
AACTGTGCGTACGGTCTGGAAACCGGCCGGTT
>d
ACTCATGTCTACGGTCCGGACCCTTACTATATAT
>e_no_barcode_match
TTTTGTCCGGACCCTTACTATATAT
>d_primer_error
AGAGTCCTGAGCGGTCCGGTACGTTTACTGGA
""".split('\n')


expected_fasta_extra_bc = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nCCCTTATATATATAT\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nCCCTTTCCA\n>Unassigned_2 c orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0\nAACCGGCCGGTT\n>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\nCCCTTACTATATAT"""

in_seqs_fixed_len_bc1_qual_scores = """>a
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 33 33 33 36 37 37 37 37 37 37 40 40 40 39 39 38
>b
35 31 31 23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37 37
>c
37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 32 32 32 36 37 35 32 32 32 32 32 32 32 32 36 37 37
>d
36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
>e_no_barcode_match
23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37
>d_primer_error
31 31 23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37""".split('\n')

# These test data have equal length barcodes, one degenerate primer
bc_map_fixed_len_bc1 = {'ACACATGTCTAC': 's1', 'AGAGTCCTGAGC': 's2',
                        'AATCGTGACTCG': 's3'}
primer_seq_lens_fixed_len_bc1 = {'ACACATGTCTAC': {'GGTCCGGA': 8},
                                 'AGAGTCCTGAGC': {'GGTCCGGA': 8},
                                 'AATCGTGACTCG': {'GGTCCGGA': 8, 'GGTCTGGA': 8}}
all_primers_fixed_len_bc1 = {'GGTCCGGA': 8, 'GGTCTGGA': 8}

# added_demultiplex barcode test data
in_seqs_added_demultiplex = """>mintberrycrunch group=coonandfriends
AATCGTGACTCGGGTCCGGACCCTTATATATATAT
>mysterion group=coonandfriends
AGAGTCCTGAGCGGTCCGGACCCTTTCCA
>mansquito group=professorchaos
AATCGTGACTCGGGTCTGGAAACCGGCCGGTT
>mintberrycrunch group=coonandfriends
AATCGTGACTTGGGTCCGGACCCTTACTATATAT
>e_no_barcode_match group=professorchaos
TTTTGTCCGGACCCTTACTATATAT
>d_primer_error group=coonandfriends
AGAGTCCTGAGCGGTCCGGTACGTTTACTGGA
""".split('\n')

# Modified barcode map to test added_demultiplex option
bc_map_added_demultiplex = {'AATCGTGACTCG,mintberry': 's1',
                            'AGAGTCCTGAGC,myst': 's2', 'AATCGTGACTCG,mansq': 's3'}

# Modified barcode map to test added_demultiplex from 'group' in label
bc_map_added_demultiplex_group = {'AATCGTGACTCG,coonandfriends': 's1',
                                  'AGAGTCCTGAGC,coonandfriends': 's2', 'AATCGTGACTCG,coonandfriends': 's3'}


expected_fasta_fixed_len_bc1 = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0
CCCTTATATATATAT
>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
CCCTTTCCA
>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0
AACCGGCCGGTT
>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1
CCCTTACTATATAT"""

# Added demultiplex options
expected_fasta_fixed_added_demultiplex = """>s1_0 mintberrycrunch orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG,mintberry bc_diffs=0\nCCCTTATATATATAT\n>s2_1 mysterion orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC,myst bc_diffs=0\nCCCTTTCCA\n>s3_2 mansquito orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG,mansq bc_diffs=0\nAACCGGCCGGTT\n>s1_3 mintberrycrunch orig_bc=AATCGTGACTTG new_bc=AATCGTGACTCG,mintberry bc_diffs=1\nCCCTTACTATATAT"""

expected_fasta_added_demultiplex_group = """>s3_0 mintberrycrunch orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG,coonandfriends bc_diffs=0\nCCCTTATATATATAT\n>s2_1 mysterion orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC,coonandfriends bc_diffs=0\nCCCTTTCCA\n>s3_2 mintberrycrunch orig_bc=AATCGTGACTTG new_bc=AATCGTGACTCG,coonandfriends bc_diffs=1\nCCCTTACTATATAT"""

# Poor quality window results in second sequence being removed
expected_fasta_fixed_len_bc1_sliding_window = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0
CCCTTATATATATAT
>s3_1 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0
AACCGGCCGGTT
>s1_2 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1
CCCTTACTATATAT"""

expected_qual_fixed_len_bc1 = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0
33 33 36 37 37 37 37 37 37 40 40 40 39 39 38
>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
36 37 37 37 37 37 37 37 37
>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0
35 32 32 32 32 32 32 32 32 36 37 37
>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1
37 37 37 37 37 37 37 37 37 37 37 37 37 37
"""

# Will be longer because primers are no longer sliced off or checked for
# mismatches
expected_fasta_fixed_len_bc1_no_primers = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0
GGTCCGGACCCTTATATATATAT
>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
GGTCCGGACCCTTTCCA
>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0
GGTCTGGAAACCGGCCGGTT
>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1
GGTCCGGACCCTTACTATATAT
>s2_4 d_primer_error orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
GGTCCGGTACGTTTACTGGA"""
# Equal length barcodes, primers, should give different results
# Due to parameter changes regarding barcode changes, primer mismatches

expected_fasta_fixed_len_bc2 = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0
CCCTTATATATATAT
>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
CCCTTTCCA
>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0
AACCGGCCGGTT
>s2_3 d_primer_error orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
ACGTTTACTGGA"""

# These test data have equal length barcodes
reverse_primers_fixed_len_bc1 = {'ACACATGTCTAC': ['CTTATAT'],
                                 'AGAGTCCTGAGC': ['GCCCTTT'],
                                 'AATCGTGACTCG': ['AGTACC']}

# Fixed barcode, reverse primers test data
in_seqs_reverse_primers = """>a
ACACATGTCTACGGTCCGGAGTACCATGATCGGCCCTTATATATATAT
>b
AGAGTCCTGAGCGGTCCGGAACCGTCGGATCAGCCCTTTCCA
>c
AATCGTGACTCGGGTCTGGAAACCGATCGACCATAGTACCGGCCGGTT
>d
ACTCATGTCTACGGTCCGGAATACGTTACGTCCCTTACTATATAT
>e_no_barcode_match
TTTTGTCCGGACCCTTACTATATAT
>d_primer_error
AGAGTCCTGAGCGGGGAGGTACGTTTACTGGA
""".split()

# expected to find the reverse rprimers for seqs a, b, and c, but does not
# find the reverse primer for d so writes out whole sequence following forward
# primer.
expected_in_seqs_reverse_primers = '>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nGTACCATGATCGGCC\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nACCGTCGGATCA\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGATCGACCAT\n>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\nATACGTTACGTCCCTTACTATATAT'

# Will truncate sequence d properly if mismatch in primer allowed
expected_in_seqs_reverse_primers_mismatch_allowed = '>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nGTACCATGATCGGCC\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nACCGTCGGATCA\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGATCGACCAT\n>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\nATACGTTACGTCC'

# Fixed barcode, reverse primers test data
in_seqs_reverse_primers_mismatch = """>a
ACACATGTCTACGGTCCGGAGTACCATGATCGGCCCCTATATATATAT
>b
AGAGTCCTGAGCGGTCCGGAACCGTCGGATCAGCCCTTTCCA
>c
AATCGTGACTCGGGTCTGGAAACCGATCGACCATAGTACCGGCCGGTT
>d
ACTCATGTCTACGGTCCGGAATACGTTACGTCCCTTACATATCCAT
>e_no_barcode_match
TTTTGTCCGGACCCTTACTATATAT
>d_primer_error
AGAGTCCTGAGCGGTCCTTTACGCCCACTGGA
""".split()

# expected to find the reverse rprimers for seqs b, and c, will not find
# seq a's reverse primer due to mismatch and will write the whole sequence
expected_in_seqs_reverse_primers_mismatch = '>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nGTACCATGATCGGCCCCTATATATATAT\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nACCGTCGGATCA\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGATCGACCAT\n>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\nATACGTTACGTCCCTTACATATCCAT'

# Will not write the d sequence, as the primer mismatches.
expected_in_seqs_reverse_primers_full_remove = '>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nGTACCATGATCGGCC\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nACCGTCGGATCA\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGATCGACCAT'

# Sample seqs.fna output file
sample_seqs_fna_file = """>a1 testing
ACACATGTCTACGGTCCGGAACGACGACGAGCGAGGGTAGC
>b2 testing more
AGAGTCCTGAGCGGTCCGGAACAGACAGGGAGAGACAGAA
>a3 testing
ACAACAGACGAGTTAGACCAA
>d4
AATCGTGACTCGGGTCTGGACAGACGAGAACGAGTTACAGACCAGA"""


expected_fasta_mad = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nCCCTTATATATATAT\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGGCCGGTT\n>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\nCCCTTACTATATAT"""


class SeqQualBadTests(TestCase):

    """Tests of the SeqQualBad class"""

    def test_init(self):
        """SeqQualBad should init OK"""
        sq = SeqQualBad('Q', None)
        self.assertEqual(sq.Name, 'Q')
        self.assertEqual(sq.F, None)
        self.assertEqual(sq.FailedIds, [])

    def test_call(self):
        """SeqQualBad should apply correct function"""
        f = lambda id_, seq, qual: len(seq) > 3
        s1 = 'aa'
        s2 = 'aaaa'
        sq = SeqQualBad('Q', f)
        self.assertEqual(sq('x', s1, [1, 2, 3]), False)
        self.assertEqual(sq('y', s2, [1, 2, 3]), True)
        self.assertEqual(sq.FailedIds, ['y'])

    def test_str(self):
        """SeqQualBad should apply correct function"""
        f = lambda id_, seq, qual: len(seq) > 3
        s1 = 'aa'
        s2 = 'aaaa'
        sq = SeqQualBad('Q', f)
        self.assertEqual(sq('x', s1, [1, 2, 3]), False)
        self.assertEqual(str(sq), 'Q\t0')
        self.assertEqual(sq('y', s2, [1, 2, 3]), True)
        self.assertEqual(str(sq), 'Q\t1')


if __name__ == '__main__':
    main()
