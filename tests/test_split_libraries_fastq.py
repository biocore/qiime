#!/usr/bin/env python
# File created on 05 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from qiime.split_libraries_fastq import (
 process_fastq_single_end_read_file,
 quality_filter_sequence,
 bad_chars_from_threshold,
 get_illumina_qual_chars,
 quality_filter_sequence)

class SplitLibrariesFastqTests(TestCase):
    """ """
    
    def setUp(self):
        self.fastq1 = fastq1.split('\n')
        self.barcode_fastq1 = barcode_fastq1.split('\n')
        self.fastq1_expected_no_qual_unassigned = \
         fastq1_expected_no_qual_unassigned
        self.fastq1_expected_default = \
         fastq1_expected_default
        self.barcode_map1 = barcode_map1
        
        
    def test_process_fastq_single_end_read_file(self):
        """process_fastq_single_end_read_file functions as expected w no qual filter
        """
        actual = process_fastq_single_end_read_file(self.fastq1,
                                       self.barcode_fastq1,
                                       self.barcode_map1,
                                       store_unassigned=True,
                                       max_bad_run_length=1000,
                                       quality_threshold='',
                                       min_per_read_length=0,
                                       rev_comp=False,
                                       rev_comp_barcode=False,
                                       barcode_in_seq=False,
                                       seq_max_N=1000,
                                       start_seq_id=0)
        actual = list(actual)
        expected = self.fastq1_expected_no_qual_unassigned
        self.assertEqual(len(actual),len(expected))
        for i in range(len(expected)):
            self.assertEqual(actual[i],expected[i])
            
    def test_process_fastq_single_end_read_file_w_defaults(self):
        """process_fastq_single_end_read_file functions as expected w default filters 
        """
        actual = process_fastq_single_end_read_file(self.fastq1,
                                       self.barcode_fastq1,
                                       self.barcode_map1,
                                       min_per_read_length=45)
        actual = list(actual)
        expected = self.fastq1_expected_default
        self.assertEqual(len(actual),len(expected))
        for i in range(len(expected)):
            self.assertEqual(actual[i],expected[i])
            
    def test_process_fastq_single_end_read_file_toggle_store_unassigned(self):
        """process_fastq_single_end_read_file handles store_unassigned
        """
        fastq_f = [
         "@990:2:4:11272:5533/1",
         "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
         "+",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
         "@990:2:4:11272:5533/2",
         "GAAAAAAAAAAT",
         "+",
         "bbbbbbbbbbbb"]
        barcode_to_sample_id = {'AAAAAAAAAAAA':'s1'}
        # empty results when store_unassigned=False
        actual = process_fastq_single_end_read_file(fastq_f,
                                                    barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    quality_threshold='B',
                                                    min_per_read_length=75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    barcode_in_seq=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)                       
        actual = list(actual)
        expected = []
        self.assertEqual(actual,expected)
        
        # non-empty results when store_unassigned=True
        actual = process_fastq_single_end_read_file(fastq_f,
                                                    barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=True,
                                                    max_bad_run_length=0,
                                                    quality_threshold='B',
                                                    min_per_read_length=75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    barcode_in_seq=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)                       
        actual = list(actual)
        expected = [('Unassigned_0 990:2:4:11272:5533/1',
         "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
         0)]
        self.assertEqual(actual,expected)
        
    def test_process_fastq_single_end_read_file_toggle_rev_comp(self):
        """process_fastq_single_end_read_file handles rev_comp
        """
        fastq_f = [\
         "@990:2:4:11272:5533/1",
         "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
         "+",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
         "@990:2:4:11272:5533/2",
         "AAAAAAAAAAAA",
         "+",
         "bbbbbbbbbbbb"]
        barcode_to_sample_id = {'AAAAAAAAAAAA':'s1'}
        
        # rev_comp = False
        actual = process_fastq_single_end_read_file(fastq_f,barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    quality_threshold='B',
                                                    min_per_read_length=75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    barcode_in_seq=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)                       
        actual = list(actual)
        expected = [('s1_0 990:2:4:11272:5533/1',
         "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
         0)]
        self.assertEqual(actual,expected)
        
        # rev_comp = True 
        actual = process_fastq_single_end_read_file(fastq_f,barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    quality_threshold='B',
                                                    min_per_read_length=75,
                                                    rev_comp=True,
                                                    rev_comp_barcode=False,
                                                    barcode_in_seq=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)                       
        actual = list(actual)
        expected = [('s1_0 990:2:4:11272:5533/1',
         "GGCTGGCTCCCCTTTCGGGGTTACCTCACCGACTTCGGGTGTTGCCGACTCTCGTGGTGTGACGGGCGGTGTGTGC",
         "`U^RY^QTTWIb_^b]aa_ab[_`a`babbbb`bbbbbbbbbbbbb`\``Ybbbbbbbbbbbbbbbbbbbbbbbbb",
         0)]
        self.assertEqual(actual,expected)
        
    def test_process_fastq_single_end_read_file_toggle_rev_comp_barcode(self):
        """process_fastq_single_end_read_file handles rev_comp_barcode
        """
        fastq_f = [\
         "@990:2:4:11272:5533/1",
         "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
         "+",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
         "@990:2:4:11272:5533/2",
         "TTTTTTTTTTTT",
         "+",
         "bbbbbbbbbbbb"]
        barcode_to_sample_id = {'AAAAAAAAAAAA':'s1'}
        # empty results when rev_comp_barcode=False
        actual = process_fastq_single_end_read_file(fastq_f,barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    quality_threshold='B',
                                                    min_per_read_length=75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    barcode_in_seq=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)                       
        actual = list(actual)
        expected = []
        self.assertEqual(actual,expected)
        
        # non-empty results when rev_comp_barcode=True
        actual = process_fastq_single_end_read_file(fastq_f,barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    quality_threshold='B',
                                                    min_per_read_length=75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=True,
                                                    barcode_in_seq=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)                       
        actual = list(actual)
        expected = [('s1_0 990:2:4:11272:5533/1',
         "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
         0)]
        self.assertEqual(actual,expected)
        
        # forward orientation no longer matches when rev_comp_barcode=True
        barcode_to_sample_id = {'TTTTTTTTTTTT':'s1'}
        actual = process_fastq_single_end_read_file(fastq_f,barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    quality_threshold='B',
                                                    min_per_read_length=75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=True,
                                                    barcode_in_seq=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = []
        self.assertEqual(actual,expected)

    def test_bad_chars_from_threshold(self):
        """bad_chars_from_threshold selects correct chars as bad 
        """
        self.assertEqual(bad_chars_from_threshold('B'),
                         {}.fromkeys(list('@AB')))
        self.assertEqual(bad_chars_from_threshold(''),{})
        self.assertEqual(bad_chars_from_threshold('~'),
                         {}.fromkeys(list(get_illumina_qual_chars())))
        self.assertEqual(bad_chars_from_threshold('@'),{'@':None})
        
    def test_quality_filter_sequence_pass(self):
        """quality_filter_sequence functions as expected for good read
        """
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))
         
    def test_quality_filter_illumina_qual(self):
        """quality_filter_sequence functions as expected with bad illumina qual digit
        """
        # header with no qual data passes
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))
         
        # header with no qual data passes
        header = "990:2:4:11271:5323/0"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))
         
        # header with no qual data passes (old barcode in header format)
        header = "HWI-6X_9267:1:1:4:1699#ACCACCC/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))
        
        # bad qual fails filter
        header = "@HWI-ST753_50:6:1101:1138:1965#0/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(3,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # bad qual passes filter if filter turned off
        header = "@HWI-ST753_50:6:1101:1138:1965#0/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=False)
        self.assertEqual(actual,(0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # good qual passes filter
        header = "@HWI-ST753_50:6:1101:1138:1965#1/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

    def test_quality_filter_sequence_fail_w_B(self):
        """quality_filter_sequence handles bad qual score as expected
        """
        
        # early 'B' in sequence causes truncation and too short of a read
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(1,"GCACTCACCGCCCGTCAC","bbbbbbbbbbbbbbbbbb"))
        
        # increasing max_bad_run_length rescues read
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=1,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))
        
        # changing threshold rescues read
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='A',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))
        
        # changing min_per_read_length rescues read
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=5,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual,(0,"GCACTCACCGCCCGTCAC","bbbbbbbbbbbbbbbbbb"))
        
    def test_quality_filter_sequence_fail_w_N(self):
        """quality_filter_sequence handles N as expected
        """
        
        # 'N' in sequence causes failure
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTNGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        expected = (2,
         "GCACTCACCGCCCGTCACACCACGAAAGTNGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`")
        self.assertEqual(actual,expected)
        
        # increasing max N rescues sequence
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTNGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=75,
                                         seq_max_N=1,
                                         filter_bad_illumina_qual_digit=True)

        expected = (0,
         "GCACTCACCGCCCGTCACACCACGAAAGTNGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`")
        self.assertEqual(actual,expected)
        
        # truncation of N rescues sequence (sequence is truncated when 
        # the quality hits B, and the truncated sequence is above the 
        # length threshold and no longer contains an N)
        header = "990:2:4:11271:5323/1"
        sequence = \
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTN"
        quality =  \
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^B`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         quality_threshold='B',
                                         min_per_read_length=50,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)

        expected = (0,
         "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTG",
         "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^")
        self.assertEqual(actual,expected)

barcode_map1 = {'AAAAAAAAAAAA':'s1',
               'AAAAAAAAAAAC':'s2',
               'AAAAAAAAAAAG':'s3',
               'AAAAAAAAAAAT':'s4',}

fastq1 = """@990:2:4:11271:5323/1
GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC
+
bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`
@990:2:4:11271:5323/1
GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG
+
bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT
@990:2:4:11272:9538/1
GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG
+
b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa
@990:2:4:11272:9538/1
GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC
+
bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O
@990:2:4:11272:7447/1
GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG
+
b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:7447/1
GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC
+
b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:19991/1
GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG
+
bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:19991/1
GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA
+
bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOcBBBBBBBBBBBBBBBBB
@990:2:4:11272:4315/1
GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG
+
bbbb_bbbbbbbbbb```Q```BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:4315/1
GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG
+
``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:5533/1
GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC
+
``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:5533/1
GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC
+
bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""

barcode_fastq1 = """@990:2:4:11271:5323/2
AAAAAAAAAAAA
+
bbbbbbbbbbbb
@990:2:4:11271:5323/2 
AAAAAAAAAAAC
+
bbcbbbbbbbbb
@990:2:4:11272:9538/2
AAAAAAAAAAAA
+
b_bbbbbbbbbb
@990:2:4:11272:9538/2
AAAAAAAAAAAT
+
bb^bbbbbbbbb
@990:2:4:11272:7447/2
AAAAAAAAAAAA
+
b`bbbbbbbbbb
@990:2:4:11272:7447/2
AAAAAAAAAAAA
+
b`bbbbbbbbbb
@990:2:4:11272:19991/2
AAAAAAAAAAAC
+
bbbbbbbbbbbb
@990:2:4:11272:19991/2
AAAAAAAAAAAC
+
bbbbbbbbbbbb
@990:2:4:11272:4315/2
AAAAAAAAAAAT
+
bbbb_bbbbbbb
@990:2:4:11272:4315/2
AAAAAAAAAAAT
+
``Q``````_``
@990:2:4:11272:5533/2 
GAAAAAAAAAAT
+
``Q``````_``
@990:2:4:11272:5533/2
AAAAAAAAAAAT
+
bbbbbbbbbbbb
"""

fastq1_expected_no_qual_unassigned = [
 ("s1_0 990:2:4:11271:5323/1",
  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
  0),
 ("s2_1 990:2:4:11271:5323/1",
  "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
  "bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT",
  1),
 ("s1_2 990:2:4:11272:9538/1",
  "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
  "b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa",
  2),
 ("s4_3 990:2:4:11272:9538/1",
  "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
  "bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O",
  3),
 ("s1_4 990:2:4:11272:7447/1",
  "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG",
  "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB",
  4),
 ("s1_5 990:2:4:11272:7447/1",
  "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC",
  "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB",
  5),
 ("s2_6 990:2:4:11272:19991/1",
  "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG",
  "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
  6),
 ("s2_7 990:2:4:11272:19991/1",
  "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA",
  "bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOcBBBBBBBBBBBBBBBBB",
  7),
 ("s4_8 990:2:4:11272:4315/1",
  "GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG",
  "bbbb_bbbbbbbbbb```Q```BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
  8),
 ("s4_9 990:2:4:11272:4315/1",
  "GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG",
  "``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
  9),
 ("Unassigned_10 990:2:4:11272:5533/1",
  "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
  "``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
  10),
 ("s4_11 990:2:4:11272:5533/1",
  "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC",
  "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
  11)]

fastq1_expected_default = [
 ("s1_0 990:2:4:11271:5323/1",
  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
  0),
 ("s2_1 990:2:4:11271:5323/1",
  "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
  "bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT",
  1),
 ("s1_2 990:2:4:11272:9538/1",
  "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
  "b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa",
  2),
 ("s4_3 990:2:4:11272:9538/1",
  "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
  "bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O",
  3),
 ("s1_4 990:2:4:11272:7447/1",
  "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
  "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`",
  4),
 ("s1_5 990:2:4:11272:7447/1",
  "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
  "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`",
  5),
 ("s2_6 990:2:4:11272:19991/1",
  "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
  "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb",
  6),
 ("s2_7 990:2:4:11272:19991/1",
  "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
  "bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc",
  7),
 ("s4_8 990:2:4:11272:5533/1",
  "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
  "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb",
  8)]


if __name__ == "__main__":
    main()