#!/usr/bin/env python
#file test_split_libraries.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Rob Knight", "William Walters"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from cogent import DNA
from StringIO import StringIO
from numpy import array
from qiime.split_libraries import (
    expand_degeneracies, get_infile, count_mismatches,
    ok_mm_primer, check_map, fasta_ids,
    count_ambig, split_seq, primer_exceeds_mismatches,
    check_barcode, make_histograms, SeqQualBad,
    seq_exceeds_homopolymers, check_window_qual_scores, check_seqs,
    local_align_primer_seq
)

class FakeOutFile(object):
    
    def __init__(self):
        self.data = ""
    
    def write(self,s):
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
         
         

    def test_check_window_qual_scores(self):
        """check_window_qual_scores returns False if window below qual 
        threshold."""
        scores1 = [8,8,8,8,8,8,8,2,2,2,2,2]
        self.assertEqual(check_window_qual_scores(scores1, 5, 5), False)
        self.assertEqual(check_window_qual_scores(scores1, 10, 5), True)
        # windowsize larger than qual score list works
        self.assertEqual(check_window_qual_scores(scores1, 100, 5), True)
        self.assertEqual(check_window_qual_scores([], 5, 1), True)
        #check each base  in its own window
        self.assertEqual(check_window_qual_scores(scores1, 1, 2), True)
        self.assertEqual(check_window_qual_scores(scores1, 1, 5), False)

    
    def test_expand_degeneracies(self):
        """generate_possibilities should make possible strings"""
        self.assertEqual(expand_degeneracies('ACG'), ['ACG'])
        self.assertEqual(expand_degeneracies('RGY'), 
            ['AGT', 'AGC', 'GGT', 'GGC'])

    def test_get_infile(self):
        """get_infile should return filehandle"""
        pass    #not practically testable, but obvious file I/O

    def test_count_mismatches(self):
        """count_mismatches should count mismatches correctly"""
        self.assertEqual(count_mismatches('GG','GG',10), 0)
        self.assertEqual(count_mismatches('GGG','AAA',10), 3)
        self.assertEqual(count_mismatches('GGG','AAA',1), 2)

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
        f.name='test.xls'
        headers, id_map, barcode_to_sample_id, warnings, errors, \
         primer_seqs_lens, all_primers = check_map(f)

        self.assertEqual(barcode_to_sample_id, {'AA':'x','AC':'y','GG':'z'})
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

    def test_fasta_ids(self):
        """fasta_ids should return list of ids in fasta files, no dups"""
        first = StringIO('>x\nACT\n>y\nAAA')
        first_copy = StringIO('>x\nACT\n>y\nAAA')
        second = StringIO('>a\nGGG\n>b\nCCC')
        self.assertEqual(fasta_ids([first, second]), set(['x','y','a','b']))
        first.seek(0) #need to reset so we can read it again
        self.assertRaises(ValueError, fasta_ids, [first,first_copy])

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
        primer = DNA.makeSequence('TAGC',Name='F5')
        seq = 'TAGC'
        # mismatch_count, hit_start
        expected = (0,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
        primer = DNA.makeSequence('TAGC',Name='F5')
        seq = 'TAGCCCCC'
        # mismatch_count, hit_start
        expected = (0,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
        primer = DNA.makeSequence('TAGC',Name='F5')
        seq = 'CCCTAGCCCCC'
        # mismatch_count, hit_start
        expected = (0,3)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)

        # different length primer
        primer = DNA.makeSequence('GTTTAGC',Name='F5')
        seq = 'GTTTAGC'
        # mismatch_count, hit_start
        expected = (0,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
        # reverse primer
        primer = DNA.makeSequence('GCTA',Name='R5')
        seq = 'TAGC'
        # mismatch_count, hit_start
        expected = (0,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
        primer = DNA.makeSequence('GCTC',Name='R5')
        seq = 'TAGCCCCC'
        # mismatch_count, hit_start
        expected = (1,2)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
        primer = DNA.makeSequence('GCTA',Name='R5')
        seq = 'CCCTAGCCCCC'
        # mismatch_count, hit_start
        expected = (1,1)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
    def test_local_align_primer_seq_fwd_rev_match_ambig(self):
        "local_align function can handle fwd/rev primers with ambigs"
        primer = DNA.makeSequence('TASC',Name='F5')
        seq = 'TAGC'
        # primer_hit, target, mismatch_count, hit_start
        expected = (0,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)

    def test_local_align_primer_seq_mm(self):
        "local_align function can handle fwd/rev primers with mismatches"
        # forward primer
        primer = DNA.makeSequence('AAAAACTTTTT',Name='F5')
        seq = 'AAAAAGTTTTT'
        # mismatch_count, hit_start
        expected = (1,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
        # forward primer
        primer = DNA.makeSequence('AAAACCTTTTT',Name='F5')
        seq = 'AAAAAGTTTTT'
        # mismatch_count, hit_start
        expected = (2,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
    def test_local_align_primer_seq_indels_middle(self):
        "local_align function can handle fwd/rev primers with indels in middle of seq"

        # Insertion in target sequence
        primer = DNA.makeSequence('CGAATCGCTATCG',Name='F5')
        seq = 'CGAATCTGCTATCG'
        # mismatch count, hit_start
        expected = (1,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        
        # Deletion in target sequence
        primer = DNA.makeSequence('CGAATCGCTATCG',Name='F5')
        seq = 'CGAATGCTATCG'
        # mismatch_count, hit_start
        expected = (1,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)

        
    def test_local_align_primer_seq_multiple_mismatch_indel(self):
        "local_align function can handle fwd/rev primers with indels and mismatches"
        # multiple insertions
        primer = DNA.makeSequence('ATCGGGCGATCATT',Name='F5')
        seq = 'ATCGGGTTCGATCATT'
        # mismatch_count, hit_start
        expected = (2,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)

        
        # two deletions
        primer = DNA.makeSequence('ACGGTACAGTGG',Name='F5')
        seq = 'ACGGCAGTGG'
        # mismatch_count, hit_start
        expected = (2,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
  
        # deletion and mismatch
        primer = DNA.makeSequence('CATCGTCGATCA',Name='F5')
        seq = 'CCTCGTGATCA'
        # mismatch_count, hit_start
        expected = (2,0)
        actual = local_align_primer_seq(primer,seq)
        self.assertEqual(actual,expected)
        

    def test_seq_exceeds_homopolymers(self):
        """seq_exceeds_homopolymers returns True if too many homopolymers"""
        self.assertEqual(seq_exceeds_homopolymers('AAACGA',3), False)
        self.assertEqual(seq_exceeds_homopolymers('AAACGA',2), True)
        self.assertEqual(seq_exceeds_homopolymers('AAACGA',1), True)
        self.assertEqual(seq_exceeds_homopolymers('AAACGATTTT',3), True)

    def test_check_barcode(self):
        """check_barcode should return False if barcode ok, True otherwise"""
        self.assertEqual(check_barcode('AA', None, ['AA']), (False, 'AA', \
         False))
        self.assertEqual(check_barcode('GCATCGTCCACA', 'golay_12', 
            ['GCATCGTCAACA']), (2, 'GCATCGTCAACA', True))
        # num errors for golay code is currently in bits
        self.assertEqual(check_barcode('AA', None, ['TT']), (True, None, False))

    def test_make_histograms(self):
        """make_histograms should make correct histograms"""
        pre_lengths = [100, 110, 105, 130, 135]
        post_lengths = [130, 135]
        pre_hist, post_hist, bin_edges = \
            make_histograms(pre_lengths, post_lengths)
        self.assertEqual(pre_hist, array([2,1,0,2]))
        self.assertEqual(post_hist, array([0,0,0,2]))
        self.assertEqual(bin_edges, array([100,110,120,130,140]))
        
    def test_check_seqs_variable_len_bc(self):
        """check_seqs handles variable length barcodes """
        
        # Simple test with variable length primers
        in_seqs = self.in_seqs_variable_len_bc1
        bc_map = self.bc_map_variable_len_bc1
        primer_seq_lens = self.primer_seq_lens_variable_len_bc1
        all_primers = self.all_primers_variable_len_bc1
        expected = self.expected_fasta_variable_len_bc1
        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=None, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="variable_length", 
         max_bc_errors=0,
         remove_unassigned=True, 
         attempt_bc_correction=False,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=0,
         disable_primer_check=False,
         reverse_primers = 'disable',
         rev_primers = {})
         
        self.assertEqual(out_f.data,expected)
        
        # Second test, includes truncated form of one of the barcodes-the
        # longest barcode should be found first
        in_seqs = self.in_seqs_variable_len_bc2
        bc_map = self.bc_map_variable_len_bc2
        primer_seq_lens = self.primer_seq_lens_variable_len_bc2
        all_primers = self.all_primers_variable_len_bc2
        expected = self.expected_fasta_variable_len_bc2
        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=None, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="variable_length", 
         max_bc_errors=0,
         remove_unassigned=True, 
         attempt_bc_correction=False,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=0,
         disable_primer_check=False,
         reverse_primers = 'disable',
         rev_primers = {})
         
        self.assertEqual(out_f.data,expected)
        

    def test_check_seqs_fixed_len_bc(self):
        """check_seqs handles fixed length barcodes """
        
        # Third test, fixed length barcodes, fixed length primers + one
        # degenerate test.  Should correct one of the passed barcodes
        in_seqs = self.in_seqs_fixed_len_bc1
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_fasta_fixed_len_bc1

        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=12, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="golay_12", 
         max_bc_errors=1.5,
         remove_unassigned=True, 
         attempt_bc_correction=True,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=0,
         disable_primer_check=False,
         reverse_primers = 'disable',
         rev_primers = {})
         
        self.assertEqual(out_f.data,expected)
        
        # Fourth test-set max_bc_errors to 0, and allow some primer mismatches
        in_seqs = self.in_seqs_fixed_len_bc2
        bc_map = self.bc_map_fixed_len_bc2
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc2
        all_primers = self.all_primers_fixed_len_bc2
        expected = self.expected_fasta_fixed_len_bc2

        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=12, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="golay_12", 
         max_bc_errors=0.5,
         remove_unassigned=True, 
         attempt_bc_correction=True,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=1,
         disable_primer_check=False,
         reverse_primers = 'disable',
         rev_primers = {})
         
        self.assertEqual(out_f.data,expected)
        

        
    def test_check_seqs_no_primers(self):
        """check_seqs handles disabled primers """
        
        # Fifth test, no primers, fixed length barcodes
        # Should correct one of the passed barcodes
        in_seqs = self.in_seqs_fixed_len_bc1
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = {}
        all_primers = {}
        expected = self.expected_fasta_fixed_len_bc1_no_primers

        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=12, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="golay_12", 
         max_bc_errors=1.5,
         remove_unassigned=True, 
         attempt_bc_correction=True,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=0,
         disable_primer_check=True,
         reverse_primers = 'disable',
         rev_primers = {})
         
        self.assertEqual(out_f.data,expected)
        
    def test_check_seqs_reverse_primers(self):
        """check_seqs handles truncating reverse primers """
        
        # Initial test, should truncate all seqs
        in_seqs = self.in_seqs_reverse_primers
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers
        rev_primers_test = self.reverse_primers_fixed_len_bc1
        
        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=12, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="golay_12", 
         max_bc_errors=1.5,
         remove_unassigned=True, 
         attempt_bc_correction=True,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=0,
         disable_primer_check=False,
         reverse_primers = 'truncate_only',
         rev_primers = rev_primers_test)
         
        self.assertEqual(out_f.data,expected)
        
        # Second test with a mismatch in seq a, should not find reverse primer
        # and will write out entire sequence.
        
        in_seqs = self.in_seqs_reverse_primers_mismatch
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers_mismatch
        rev_primers_test = self.reverse_primers_fixed_len_bc1
        
        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=12, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="golay_12", 
         max_bc_errors=1.5,
         remove_unassigned=True, 
         attempt_bc_correction=True,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=0,
         disable_primer_check=False,
         reverse_primers = 'truncate_only',
         rev_primers = rev_primers_test)
         
        self.assertEqual(out_f.data,expected)
        
        # With mismatches allowed set to 1, should restore truncation.
        in_seqs = self.in_seqs_reverse_primers_mismatch
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers_mismatch_allowed
        rev_primers_test = self.reverse_primers_fixed_len_bc1
        
        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=12, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="golay_12", 
         max_bc_errors=1.5,
         remove_unassigned=True, 
         attempt_bc_correction=True,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=1,
         disable_primer_check=False,
         reverse_primers = 'truncate_only',
         rev_primers = rev_primers_test)
         
        self.assertEqual(out_f.data,expected)
        
        # Testing truncate_remove, which should not write sequences where
        # the reverse primer is not found
        in_seqs = self.in_seqs_reverse_primers
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers_full_remove
        rev_primers_test = self.reverse_primers_fixed_len_bc1
        

        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=12, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="golay_12", 
         max_bc_errors=1.5,
         remove_unassigned=True, 
         attempt_bc_correction=True,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=0,
         disable_primer_check=False,
         reverse_primers = 'truncate_remove',
         rev_primers = rev_primers_test)
         
        self.assertEqual(out_f.data,expected)
        
        # Testing truncate_remove, with mismatch set to 1 should allow
        # all 4 sequences to be written, truncated
        in_seqs = self.in_seqs_reverse_primers_mismatch
        bc_map = self.bc_map_fixed_len_bc1
        primer_seq_lens = self.primer_seq_lens_fixed_len_bc1
        all_primers = self.all_primers_fixed_len_bc1
        expected = self.expected_in_seqs_reverse_primers_mismatch_allowed
        rev_primers_test = self.reverse_primers_fixed_len_bc1
        

        
        
        out_f = FakeOutFile()
        
        actual = check_seqs(
         fasta_out=out_f, 
         fasta_files = [in_seqs], 
         starting_ix=0, 
         valid_map = bc_map, 
         qual_mappings={}, 
         filters=[], 
         barcode_len=12, 
         keep_primer=False, 
         keep_barcode=False, 
         barcode_type="golay_12", 
         max_bc_errors=1.5,
         remove_unassigned=True, 
         attempt_bc_correction=True,
         primer_seqs_lens=primer_seq_lens,
         all_primers=all_primers, 
         max_primer_mm=1,
         disable_primer_check=False,
         reverse_primers = 'truncate_remove',
         rev_primers = rev_primers_test)
         
        self.assertEqual(out_f.data,expected)
        

        



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

bc_map_variable_len_bc1 = {'ACC':'s1','AGGA':'s2','ATTA':'s3'}
primer_seq_lens_variable_len_bc1 = {'ACC':{'GGTCCGGA':8},
                                    'AGGA':{'GTCCGGA':7},
                                    'ATTA':{'ACCCGGA':7}}
all_primers_variable_len_bc1 = {'GGTCCGGA':8,'GTCCGGA':7,'ACCCGGA':7}

expected_fasta_variable_len_bc1 = """>s1_0 a orig_bc=ACC new_bc=ACC bc_diffs=0
CCCTTATATATATAT
>s2_1 b orig_bc=AGGA new_bc=AGGA bc_diffs=0
CCCTTTCCA
>s3_2 c orig_bc=ATTA new_bc=ATTA bc_diffs=0
AACCGGCCGGTT
>s1_3 d orig_bc=ACC new_bc=ACC bc_diffs=0
CCCTTACTATATAT
"""
bc_map_variable_len_bc2 = {'ACC':'s1','AGGA':'s2','ATTA':'s3','AGG':'s4'}
primer_seq_lens_variable_len_bc2 = {'ACC':{'GGTCCGGA':8},
                                    'AGGA':{'GTCCGGA':7},
                                    'ATTA':{'ACCCGGA':7},
                                    'AGG':{'AGTCCGGA':8}}
all_primers_variable_len_bc2 = {'GGTCCGGA':8,'GTCCGGA':7,'ACCCGGA':7,
 'AGTCCGGA':8}
 
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
""".split()

# These test data have equal length barcodes, one degenerate primer
bc_map_fixed_len_bc1 = {'ACACATGTCTAC':'s1','AGAGTCCTGAGC':'s2',
 'AATCGTGACTCG':'s3'}
primer_seq_lens_fixed_len_bc1 = {'ACACATGTCTAC':{'GGTCCGGA':8},
                                    'AGAGTCCTGAGC':{'GGTCCGGA':8},
                                    'AATCGTGACTCG':{'GGTCCGGA':8,'GGTCTGGA':8}}
all_primers_fixed_len_bc1 = {'GGTCCGGA':8,'GGTCTGGA':8}


expected_fasta_fixed_len_bc1 = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0
CCCTTATATATATAT
>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
CCCTTTCCA
>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0
AACCGGCCGGTT
>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1
CCCTTACTATATAT
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
GGTCCGGTACGTTTACTGGA
"""
# Equal length barcodes, primers, should give different results
# Due to parameter changes regarding barcode changes, primer mismatches

expected_fasta_fixed_len_bc2 = """>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0
CCCTTATATATATAT
>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
CCCTTTCCA
>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0
AACCGGCCGGTT
>s2_3 d_primer_error orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0
ACGTTTACTGGA
"""

# These test data have equal length barcodes
reverse_primers_fixed_len_bc1 = {'ACACATGTCTAC':'CTTATAT',
                                    'AGAGTCCTGAGC':'GCCCTTT',
                                    'AATCGTGACTCG':'AGTACC'}

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
expected_in_seqs_reverse_primers = '>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nGTACCATGATCGGCC\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nACCGTCGGATCA\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGATCGACCAT\n>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\nATACGTTACGTCCCTTACTATATAT\n'

# Will truncate sequence d properly if mismatch in primer allowed
expected_in_seqs_reverse_primers_mismatch_allowed = '>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nGTACCATGATCGGCC\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nACCGTCGGATCA\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGATCGACCAT\n>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\nATACGTTACGTCC\n'

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
expected_in_seqs_reverse_primers_mismatch = '>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nGTACCATGATCGGCCCCTATATATATAT\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nACCGTCGGATCA\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGATCGACCAT\n>s1_3 d orig_bc=ACTCATGTCTAC new_bc=ACACATGTCTAC bc_diffs=1\nATACGTTACGTCCCTTACATATCCAT\n'

# Will not write the d sequence, as the primer mismatches.
expected_in_seqs_reverse_primers_full_remove = '>s1_0 a orig_bc=ACACATGTCTAC new_bc=ACACATGTCTAC bc_diffs=0\nGTACCATGATCGGCC\n>s2_1 b orig_bc=AGAGTCCTGAGC new_bc=AGAGTCCTGAGC bc_diffs=0\nACCGTCGGATCA\n>s3_2 c orig_bc=AATCGTGACTCG new_bc=AATCGTGACTCG bc_diffs=0\nAACCGATCGACCAT\n'


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
        self.assertEqual(sq('x',s1, [1,2,3]), False)
        self.assertEqual(sq('y',s2, [1,2,3]), True)
        self.assertEqual(sq.FailedIds, ['y'])

    def test_str(self):
        """SeqQualBad should apply correct function"""
        f = lambda id_, seq, qual: len(seq) > 3
        s1 = 'aa'
        s2 = 'aaaa'
        sq = SeqQualBad('Q', f)
        self.assertEqual(sq('x',s1, [1,2,3]), False)
        self.assertEqual(str(sq), 'Q\t0')
        self.assertEqual(sq('y',s2, [1,2,3]), True)
        self.assertEqual(str(sq), 'Q\t1')
        



if __name__ =='__main__':
    main()
