#!/usr/bin/env python
# File created on 22 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from cogent.util.unit_test import TestCase, main
from qiime.util import get_tmp_filename
from cogent.util.misc import remove_files
from qiime.parse import parse_mapping_file, IlluminaParseError
from qiime.split_libraries_illumina import (
    parse_illumina_paired_end_read_files,
    parse_illumina_single_end_read_file,
    mapping_data_to_barcode_map,
    read_qual_score_filter, 
    bad_chars_from_threshold,
    process_illumina_paired_end_read_files,
    illumina_read_description_from_read_data,
    process_illumina_single_end_read_file)

class IlluminaParserTests(TestCase):
    """
    """
    
    def setUp(self):
        """ """
        
        self.files_to_remove = []
        
        self.illumina_read1 = illumina_read1
        self.illumina_read2 = illumina_read2
        self.illumina_read1_seq_N = illumina_read1_seq_N
        self.illumina_read2_seq_N = illumina_read2_seq_N
        
        self.mapping_f = mapping_f

        self.expected_seqs_file1 = expected_seqs_file1
        self.expected_qual_file1 = expected_qual_file1
        self.barcode_to_sample_id1 = barcode_to_sample_id1
        
        self.expected_seqs_file2 = expected_seqs_file2
        self.expected_qual_file2 = expected_qual_file2
        self.barcode_to_sample_id2 = barcode_to_sample_id2
        
        self.expected_seqs_file3 = expected_seqs_file3
        self.expected_qual_file3 = expected_qual_file3
        self.barcode_to_sample_id3 = barcode_to_sample_id3
        
        self.expected_5prime_seqs_file1 = expected_5prime_seqs_file1
        self.expected_5prime_qual_file1 = expected_5prime_qual_file1
        self.expected_5prime_seqs_file2 = expected_5prime_seqs_file2
        self.expected_5prime_qual_file2 = expected_5prime_qual_file2
        
        self.barcode_to_sample_id4 = barcode_to_sample_id4
        self.illumina_read3_bc_in_seq = illumina_read3_bc_in_seq
        self.expected_seq_file_bc_in_seq1 = expected_seq_file_bc_in_seq1
        self.expected_qual_file_bc_in_seq1 = expected_qual_file_bc_in_seq1
    
    def tearDown(self):
        remove_files(self.files_to_remove)
    
    def test_illumina_read_description_from_read_data(self):
        """illumina_read_description_from_read_data: functions as expected
        """
        d = {\
         'Machine Name':'HWI-6X_9267',\
         'Channel Number':1,\
         'Tile Number':1,\
         'X Position':4,\
         'Y Position':390,\
         'Barcode':'GGAGGT',\
         'Full Y Position Field':'390#ACCTCCC/1',\
         'Sequence':\
          'GGTGGNTGATTGGANTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',\
         'Quality Score':\
          'BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBAAAAAAAAAAABBBB'}
        actual = illumina_read_description_from_read_data(d)
        expected = 'HWI-6X_9267:1:1:4:390#ACCTCCC'
        self.assertEqual(actual,expected)
        
    def test_parse_illumina_single_end_read_file_no_revComp(self):
        """parse_illumina_single_end_read_file: single end read parsing functions as expected (+ strand)
        """
        actual = list(parse_illumina_single_end_read_file(illumina_read1,barcode_length=6,\
            max_bad_run_length=0,quality_threshold=1e-5,min_per_read_length=70,
            rev_comp=False,rev_comp_barcode=True,barcode_in_seq=False))
        expected =[
         ('HWI-6X_9267:1:1:4:1699#ACCACCC','GGTGGT',\
          'TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCC'+\
          'CCCCCCCCCCCCCCCCCCCCCCCC',\
          'abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaa'+\
          'aaaaaaaaaaaaaaaaaaaaaaaa'),\
        ('HWI-6X_9267:1:1:4:390#ACCTCCC','GGAGGT',\
          'GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCC'+\
          'CCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAA',\
          'aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'+\
          'aaaaaaaaaaaaaaaaaaaaaaaaaaaa')]
        self.assertEqual(actual,expected)
        
    def test_parse_illumina_single_end_read_file_revComp(self):
        """parse_illumina_single_end_read_file: single end read parsing functions as expected (- strand)
        """
        actual = list(parse_illumina_single_end_read_file(illumina_read2,barcode_length=6,\
            max_bad_run_length=0,quality_threshold=1e-5,min_per_read_length=70,
            rev_comp=True,rev_comp_barcode=True,barcode_in_seq=False))
        expected =[
         ('HWI-6X_9267:1:1:4:1699#ACCACCC','GGTGGT',\
          'GGTTTTTTTTTAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG'+\
          'GGGGCCCCCCCCCCCTTTTTTTTAAAA',\
          'cccccccccccccccccbbbbbbbbbbbbbbbbbbbbbbbbba'+\
          'aaaaaaaaaaaaaaaaaaaaaaaaaaaa'),\
        ('HWI-6X_9267:1:1:4:390#ACCTCCC','GGAGGT',\
          'CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA'+\
          'CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',\
          'bbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbb'+\
          'bbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaa')]
        self.assertEqual(actual,expected)
     
    def test_parse_illumina_single_end_read_file_N_revComp(self):
        """parse_illumina_single_end_read_file: single end read parsing functions with N (- strand)
        """
        actual = list(parse_illumina_single_end_read_file(illumina_read2_N,barcode_length=6,\
            max_bad_run_length=0,quality_threshold=1e-5,min_per_read_length=70,
            rev_comp=True,rev_comp_barcode=True,barcode_in_seq=False))
        expected =[]
        self.assertEqual(actual,expected)
        
    def test_parse_illumina_single_end_read_file_N_revComp(self):
        """parse_illumina_single_end_read_file: single end read parsing functions with N (- strand)
        """
        actual = list(parse_illumina_single_end_read_file(illumina_read2_N,barcode_length=6,\
            max_bad_run_length=0,quality_threshold=1e-5,min_per_read_length=70,
            rev_comp=True,rev_comp_barcode=True,barcode_in_seq=False))
        expected =[]
        self.assertEqual(actual,expected)
        
        # allow one N in barcode
        actual = list(parse_illumina_single_end_read_file(illumina_read2_N,barcode_length=6,\
            max_bad_run_length=0,quality_threshold=1e-5,min_per_read_length=70,
            rev_comp=True,rev_comp_barcode=True,barcode_in_seq=False,
            barcode_max_N=1,seq_max_N=0))
        expected =[
         ('HWI-6X_9267:1:1:4:390#ACNTCCC','GGANGT',\
          'CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA'+\
          'CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',\
          'bbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbb'+\
          'bbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaa')]
        self.assertEqual(actual,expected)
        
        # allow one N in barcode
        actual = list(parse_illumina_single_end_read_file(illumina_read2_N,barcode_length=6,\
            max_bad_run_length=0,quality_threshold=1e-5,min_per_read_length=70,
            rev_comp=True,rev_comp_barcode=True,barcode_in_seq=False,
            barcode_max_N=0,seq_max_N=1))
        expected =[
         ('HWI-6X_9267:1:1:4:1699#ACCACCC','GGTGGT',\
          'GGTTTTTTTTTAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGCC'+\
          'CCCNCCCCCTTTTTTTTAAAA',\
          'cccccccccccccccccbbbbbbbbbbbbbbbbbbbbbbbbba'+\
          'aaaaaaaaaaaaaaaaaaaaaaaaaaaa')]
        self.assertEqual(actual,expected)
        
        # allow one N in barcode and one N in seq
        actual = list(parse_illumina_single_end_read_file(illumina_read2_N,barcode_length=6,\
            max_bad_run_length=0,quality_threshold=1e-5,min_per_read_length=70,
            rev_comp=True,rev_comp_barcode=True,barcode_in_seq=False,
            barcode_max_N=1,seq_max_N=1))
        expected =[
         ('HWI-6X_9267:1:1:4:1699#ACCACCC','GGTGGT',\
          'GGTTTTTTTTTAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGCC'+\
          'CCCNCCCCCTTTTTTTTAAAA',\
          'cccccccccccccccccbbbbbbbbbbbbbbbbbbbbbbbbba'+\
          'aaaaaaaaaaaaaaaaaaaaaaaaaaaa'),\
         ('HWI-6X_9267:1:1:4:390#ACNTCCC','GGANGT',\
          'CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA'+\
          'CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',\
          'bbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbb'+\
          'bbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaa')]
        self.assertEqual(actual,expected)

    def test_parse_illumina_paired_end_read_files(self):
        """parse_illumina_paired_end_read_files: functions as expected """
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1,illumina_read2,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=70))
        expected = [
         ('HWI-6X_9267:1:1:4:1699#ACCACCC','GGTGGT',\
          'TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'+\
          'CCCCCCCCCCCCCCCCCGGTTTTTTTTTAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGG'+\
          'GGGGGGCCCCCCCCCCCTTTTTTTTAAAA',\
          'abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaa'+\
          'aaaaaaaaaaaaaaaaaaaaaaaacccccccccccccccccbbbbbbbbbbbbbbbbbbbbbbb'+\
          'bbaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'),\
         ('HWI-6X_9267:1:1:4:390#ACCTCCC','GGAGGT',\
          'GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCC'+\
          'CCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAACGTACGTACGT'+\
          'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC'+\
          'GTACGTACGTACGTACGTACGTACGT',\
          'aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaa'+\
          'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbb'+\
          'bbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbb'+\
          'bbbbaaaaaaaaaaaaaaaaaaaaaaaaaa')]
        self.assertEqual(actual,expected)
    
        # alt min_per_read_length
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1,illumina_read2,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=75))
        expected = [
         ('HWI-6X_9267:1:1:4:390#ACCTCCC','GGAGGT',\
          'GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCC'+\
          'CCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAACGTACGTACGT'+\
          'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC'+\
          'GTACGTACGTACGTACGTACGTACGT',\
          'aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaa'+\
          'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbb'+\
          'bbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbb'+\
          'bbbbaaaaaaaaaaaaaaaaaaaaaaaaaa')]
        self.assertEqual(actual,expected)
    
        # alt quality_threshold (just checking number of results
        # to be sure that alt value is passed through)
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1,illumina_read2,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-55,
         min_per_read_length=75))
        self.assertEqual(len(actual),0)
        
        # alt max_bad_run_length (just checking number of results
        # to be sure that alt value is passed through)
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1,illumina_read2,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-55,
         min_per_read_length=75))
        self.assertEqual(len(actual),0)
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1,illumina_read2,barcode_length=6,\
         max_bad_run_length=150,rev_comp_barcode=True,quality_threshold=1e-55,
         min_per_read_length=75))
        self.assertEqual(len(actual),2)
        
    def test_parse_illumina_paired_end_read_files_N(self):
        """parse_illumina_paired_end_read_files: functions as expected with N chars"""
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1_N,illumina_read2,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=70))
        expected = []
        self.assertEqual(actual,expected)
        
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1,illumina_read2_N,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=70))
        expected = []
        self.assertEqual(actual,expected)
        
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1_N,illumina_read2_N,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=70))
        expected = []
        self.assertEqual(actual,expected)
        
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1_N,illumina_read2_N,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=70,barcode_max_N=1))
        self.assertEqual(len(actual),1)
        
        # one sequence discarded due to barcode N, other due to 2 sequence Ns
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1_N,illumina_read2_N,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=70,seq_max_N=1))
        self.assertEqual(len(actual),0)
        
        actual = list(parse_illumina_paired_end_read_files(\
         illumina_read1_N,illumina_read2_N,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=70,barcode_max_N=1,seq_max_N=2))
        self.assertEqual(len(actual),2)
        
    def test_parse_illumina_paired_end_read_files_error(self):
        """parse_illumina_paired_end_read_files: detects mis-matched lines """
        reversed_illumina_read_lines1 = \
         [self.illumina_read1[1],self.illumina_read1[0]]
        record_iter = parse_illumina_paired_end_read_files(\
         reversed_illumina_read_lines1,illumina_read2,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=75)
        self.assertRaises(IlluminaParseError,list,record_iter)
        
        reversed_illumina_read_lines2 = \
         [self.illumina_read2[1],self.illumina_read2[0]]
        record_iter = parse_illumina_paired_end_read_files(\
         illumina_read1,reversed_illumina_read_lines2,barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=75)
        self.assertRaises(IlluminaParseError,list,record_iter)
        
        # no error if reversed-order files are passed
        list(parse_illumina_paired_end_read_files(\
         reversed_illumina_read_lines1,reversed_illumina_read_lines2,\
         barcode_length=6,\
         max_bad_run_length=0,rev_comp_barcode=True,quality_threshold=1e-5,
         min_per_read_length=75))
    
    def test_read_qual_score_filter(self):
        """read_qual_score_filter: functions as expected
        """
        actual = read_qual_score_filter('ACGTACGTA','aa``ba^``',\
            max_run_length=0,threshold=1e-5)
        expected = ('ACGTACGTA','aa``ba^``')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','Ba``ba^``',\
            max_run_length=0,threshold=1e-5)
        expected = ('','')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','aB``ba^``',\
            max_run_length=0,threshold=1e-5)
        expected = ('A','a')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','aa``ba^`B',\
            max_run_length=0,threshold=1e-5)
        expected = ('ACGTACGT','aa``ba^`')
        
        self.assertEqual(actual,expected)
        actual = read_qual_score_filter(\
         '','',max_run_length=0,threshold=1e-5)
        expected = ('','')
        
        self.assertEqual(actual,expected)
        actual = read_qual_score_filter(\
            'ACGTACGTANNNNNNNNNAAAAAAAACCCCCCCC',\
            'aa``ba^``bbbbbbbbbbbbbbbbbbbbbbbbb',\
            max_run_length=0,threshold=1e-5)
        expected = ('ACGTACGTANNNNNNNNNAAAAAAAACCCCCCCC',\
            'aa``ba^``bbbbbbbbbbbbbbbbbbbbbbbbb')
        
    def test_read_qual_score_filter_alt_p_threshold(self):
        """read_qual_score_filter: functions as expected w alt p threshold
        """
        actual = read_qual_score_filter('ACGTACGTA','aa``ba^``',\
            max_run_length=0,threshold=1e-5)
        expected = ('ACGTACGTA','aa``ba^``')
        self.assertEqual(actual,expected)
        
        # change threshold to filter/not filter on C
        actual = read_qual_score_filter('ACGTACGTA','aa``Ca^``',\
            max_run_length=0,threshold=1e-5)
        expected = ('ACGT','aa``')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','aa``Ca^``',\
            max_run_length=0,threshold=1e-2)
        expected = ('ACGTACGTA','aa``Ca^``')
        self.assertEqual(actual,expected)
        
        # Lowest prob value is 1e-62 associated with ~
        actual = read_qual_score_filter('ACGTACGTA','~~~~~~~~~',\
            max_run_length=0,threshold=1e-62)
        expected = ('ACGTACGTA','~~~~~~~~~')
        self.assertEqual(actual,expected)
        actual = read_qual_score_filter('ACGTACGTA','~~~~~~~~~',\
            max_run_length=0,threshold=1e-63)
        expected = ('','')
        self.assertEqual(actual,expected)
        
        # Highest prob value is 1. associated with @
        actual = read_qual_score_filter('ACGTACGTA','@@@@@@@@@',\
            max_run_length=0,threshold=1.)
        expected = ('ACGTACGTA','@@@@@@@@@')
        actual = read_qual_score_filter('ACGTACGTA','@@@@@@@@@',\
            max_run_length=0,threshold=0.1)
        expected = ('','')
        
        
    def test_read_qual_score_filter_alt_bad_run_length(self):
        """read_qual_score_filter: functions as expected w alt max_run_length
        """
        actual = read_qual_score_filter('ACGTACGTA','aa``Ba^``',\
            max_run_length=0,threshold=1e-5)
        expected = ('ACGT','aa``')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','BBBBBBBBB',\
            max_run_length=0,threshold=1e-5)
        expected = ('','')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','Ba``Ba^``',\
            max_run_length=0,threshold=1e-5)
        expected = ('','')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','aa``Ba^``',\
            max_run_length=1,threshold=1e-5)
        expected = ('ACGTACGTA','aa``Ba^``')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','aa``B^B``',\
            max_run_length=1,threshold=1e-5)
        expected = ('ACGTACGTA','aa``B^B``')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','aa``BB^``',\
            max_run_length=1,threshold=1e-5)
        expected = ('ACGT','aa``')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','aa``BB^``',\
            max_run_length=2,threshold=1e-5)
        expected = ('ACGTACGTA','aa``BB^``')
        self.assertEqual(actual,expected)
        
        actual = read_qual_score_filter('ACGTACGTA','BBBBBBBBB',\
            max_run_length=9,threshold=1e-5)
        expected = ('ACGTACGTA','BBBBBBBBB')
        self.assertEqual(actual,expected)
        
    def test_bad_chars_from_threshold(self):
        """bad_chars_from_threshold: functions as expected """
        actual = bad_chars_from_threshold(1e-5)
        expected = {}.fromkeys(list('@ABCD'))
        self.assertEqual(actual,expected)
        
        actual = bad_chars_from_threshold(1e-7)
        expected = {}.fromkeys(list('@ABCDEF'))
        self.assertEqual(actual,expected)
        
        actual = bad_chars_from_threshold(1e-63)
        expected = {}.fromkeys(list('@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcd'+\
         'efghijklmnopqrstuvwxyz{|}~'))
        self.assertEqual(actual,expected)
        
        actual = bad_chars_from_threshold(1e-62)
        expected = {}.fromkeys(list('@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcd'+\
         'efghijklmnopqrstuvwxyz{|}'))
        self.assertEqual(actual,expected)
        
        actual = bad_chars_from_threshold(1.)
        expected = {}.fromkeys(list(''))
        self.assertEqual(actual,expected)
        
        
    def test_mapping_data_to_barcode_map(self):
        """parse_barcode_map: functions as expected
        """
        mapping_data, mapping_headers, mapping_comments =\
            parse_mapping_file(self.mapping_f)
        expected = {'GGTGGT':'Samp2',\
                    'GGAGGT':'SAMP_1',\
                    'GGTTAA':'dflsdflsdfsdfsdfsd'}
        self.assertEqual(mapping_data_to_barcode_map(mapping_data),expected)
        
    def test_process_illumina_single_end_read_file1(self):
        """process_illumina_single_end_read_file: functions as expected
        """
        output_seqs_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.fasta')
        output_qual_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        
        open(read_fp,'w').write('\n'.join(self.illumina_read1))
        self.files_to_remove.append(read_fp)
        
        actual = process_illumina_single_end_read_file(read_fp,output_seqs_fp,
         output_qual_fp,barcode_to_sample_id=self.barcode_to_sample_id1, 
         barcode_length=6,store_unassigned=True,max_bad_run_length=0,
         quality_threshold=1e-5,min_per_read_length=70,rev_comp=False,
         rev_comp_barcode=True,barcode_in_seq=False,seq_max_N=0,start_seq_id=0)
        
        self.files_to_remove.append(output_seqs_fp)
        self.files_to_remove.append(output_qual_fp)

        # next_seq_id is returned correctly
        self.assertEqual(actual,2)
        
        # correct seq file is returned
        self.assertEqual([l.strip() for l in list(open(output_seqs_fp))],\
         self.expected_5prime_seqs_file1)
        
        # correct qual file is returned
        self.assertEqual([l.strip() for l in list(open(output_qual_fp))],\
         self.expected_5prime_qual_file1)

        
    def test_process_illumina_single_end_read_file2(self):
        """process_illumina_single_end_read_file: alt seq max N
        """
        output_seqs_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.fasta')
        output_qual_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        
        open(read_fp,'w').write('\n'.join(self.illumina_read1_seq_N))
        self.files_to_remove.append(read_fp)
        
        ## seq_max_N=1 allows both sequences
        actual = process_illumina_single_end_read_file(read_fp,output_seqs_fp,
         output_qual_fp,barcode_to_sample_id=self.barcode_to_sample_id1, 
         barcode_length=6,store_unassigned=True,max_bad_run_length=0,
         quality_threshold=1e-5,min_per_read_length=70,rev_comp=False,
         rev_comp_barcode=True,barcode_in_seq=False,seq_max_N=1,start_seq_id=0)
        
        self.files_to_remove.append(output_seqs_fp)
        self.files_to_remove.append(output_qual_fp)

        # next_seq_id is returned correctly
        self.assertEqual(actual,2)
        
        # correct seq file is returned
        self.assertEqual([l.strip() for l in list(open(output_seqs_fp))],\
         self.expected_5prime_seqs_file2)
        
        # correct qual file is returned
        self.assertEqual([l.strip() for l in list(open(output_qual_fp))],\
         self.expected_5prime_qual_file2)
         
        ## Lower seq_max_N yields no results
        actual = process_illumina_single_end_read_file(read_fp,output_seqs_fp,
         output_qual_fp,barcode_to_sample_id=self.barcode_to_sample_id1, 
         barcode_length=6,store_unassigned=True,max_bad_run_length=0,
         quality_threshold=1e-5,min_per_read_length=70,rev_comp=False,
         rev_comp_barcode=True,barcode_in_seq=False,seq_max_N=0,start_seq_id=0)

        # next_seq_id is returned correctly
        self.assertEqual(actual,0)
        
        # correct seq file is returned
        self.assertEqual([l.strip() for l in list(open(output_seqs_fp))],\
         [])
        
        # correct qual file is returned
        self.assertEqual([l.strip() for l in list(open(output_qual_fp))],\
         [])


    def test_process_illumina_paired_end_read_files1(self):
        """process_illumina_paired_end_read_files: functions as expected
        """
        output_seqs_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.fasta')
        output_qual_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read1_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read2_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        
        open(read1_fp,'w').write('\n'.join(self.illumina_read1))
        self.files_to_remove.append(read1_fp)
        open(read2_fp,'w').write('\n'.join(self.illumina_read2))
        self.files_to_remove.append(read2_fp)
        
        actual = process_illumina_paired_end_read_files(\
         read1_fp,read2_fp,output_seqs_fp,output_qual_fp,\
         barcode_to_sample_id=self.barcode_to_sample_id1,\
         barcode_length=6,rev_comp_barcode=True,\
         store_unassigned=True,max_bad_run_length=0,\
         quality_threshold=1e-5,min_per_read_length=70,\
         start_seq_id=0)
        
        self.files_to_remove.append(output_seqs_fp)
        self.files_to_remove.append(output_qual_fp)
        
        # next_seq_id is returned correctly
        self.assertEqual(actual,2)
        
        # correct seq file is returned
        self.assertEqual([l.strip() for l in list(open(output_seqs_fp))],\
         self.expected_seqs_file1)
        
        # correct qual file is returned
        self.assertEqual([l.strip() for l in list(open(output_qual_fp))],\
         self.expected_qual_file1)
         

    def test_process_illumina_paired_end_read_files2(self):
        """process_illumina_paired_end_read_files: functions as expected with alt start_seq_id
        """
        output_seqs_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.fasta')
        output_qual_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read1_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read2_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        
        open(read1_fp,'w').write('\n'.join(self.illumina_read1))
        self.files_to_remove.append(read1_fp)
        open(read2_fp,'w').write('\n'.join(self.illumina_read2))
        self.files_to_remove.append(read2_fp)
        
        actual = process_illumina_paired_end_read_files(\
         read1_fp,read2_fp,output_seqs_fp,output_qual_fp,\
         barcode_to_sample_id=self.barcode_to_sample_id2,\
         barcode_length=6,rev_comp_barcode=True,\
         store_unassigned=True,max_bad_run_length=0,\
         quality_threshold=1e-5,min_per_read_length=70,\
         start_seq_id=42)
        
        self.files_to_remove.append(output_seqs_fp)
        self.files_to_remove.append(output_qual_fp)
        
        # next_seq_id is returned correctly
        self.assertEqual(actual,44)
        
        # correct seq file is returned
        self.assertEqual([l.strip() for l in list(open(output_seqs_fp))],\
         self.expected_seqs_file2)
        
        # correct qual file is returned
        self.assertEqual([l.strip() for l in list(open(output_qual_fp))],\
         self.expected_qual_file2)

    def test_process_illumina_paired_end_read_files3(self):
        """process_illumina_paired_end_read_files: functions as expected with seq_max_N
        """
        output_seqs_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.fasta')
        output_qual_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read1_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read2_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        
        open(read1_fp,'w').write('\n'.join(self.illumina_read1_seq_N))
        self.files_to_remove.append(read1_fp)
        open(read2_fp,'w').write('\n'.join(self.illumina_read2_seq_N))
        self.files_to_remove.append(read2_fp)
        
        # seq_max_N=2 allows both sequences
        actual = process_illumina_paired_end_read_files(\
         read1_fp,read2_fp,output_seqs_fp,output_qual_fp,\
         barcode_to_sample_id=self.barcode_to_sample_id3,\
         barcode_length=6,rev_comp_barcode=True,\
         store_unassigned=True,max_bad_run_length=0,\
         quality_threshold=1e-5,min_per_read_length=70,\
         start_seq_id=0,seq_max_N=2)
        
        self.files_to_remove.append(output_seqs_fp)
        self.files_to_remove.append(output_qual_fp)
        
        # next_seq_id is returned correctly
        self.assertEqual(actual,2)
        
        # correct seq file is returned
        self.assertEqual([l.strip() for l in list(open(output_seqs_fp))],\
         self.expected_seqs_file3)
        
        # correct qual file is returned
        self.assertEqual([l.strip() for l in list(open(output_qual_fp))],\
         self.expected_qual_file3)

        # Lower seq_max_N returns no results
        actual = process_illumina_paired_end_read_files(\
         read1_fp,read2_fp,output_seqs_fp,output_qual_fp,\
         barcode_to_sample_id=self.barcode_to_sample_id3,\
         barcode_length=6,rev_comp_barcode=True,\
         store_unassigned=True,max_bad_run_length=0,\
         quality_threshold=1e-5,min_per_read_length=70,\
         start_seq_id=0,seq_max_N=0)
        
        # next_seq_id is returned correctly
        self.assertEqual(actual,0)
        
        # correct seq file is returned
        self.assertEqual([l.strip() for l in list(open(output_seqs_fp))],[])
        
        # correct qual file is returned
        self.assertEqual([l.strip() for l in list(open(output_qual_fp))],[])
        
    def test_process_illumina_single_end_read_file_bc_in_seq(self):
        """process_illumina_single_end_read_file: functions as expected with bc in seq
        """
        output_seqs_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.fasta')
        output_qual_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        read_fp = get_tmp_filename(\
         prefix='ParseIlluminaTests',suffix='.txt')
        
        open(read_fp,'w').write('\n'.join(self.illumina_read3_bc_in_seq))
        self.files_to_remove.append(read_fp)
        
        actual = process_illumina_single_end_read_file(read_fp,output_seqs_fp,
         output_qual_fp,barcode_to_sample_id=self.barcode_to_sample_id4, 
         barcode_length=12,store_unassigned=True,max_bad_run_length=0,
         quality_threshold=1e-5,min_per_read_length=70,rev_comp=False,
         rev_comp_barcode=False,barcode_in_seq=True,
         seq_max_N=0,start_seq_id=0)
        
        self.files_to_remove.append(output_seqs_fp)
        self.files_to_remove.append(output_qual_fp)

        # next_seq_id is returned correctly
        self.assertEqual(actual,2)
        
        # correct seq file is returned
        self.assertEqual([l.strip() for l in list(open(output_seqs_fp))],\
         self.expected_seq_file_bc_in_seq1)
        
        # correct qual file is returned
        self.assertEqual([l.strip() for l in list(open(output_qual_fp))],\
         self.expected_qual_file_bc_in_seq1)

illumina_read1 = """HWI-6X_9267:1:1:4:1699#ACCACCC/1:TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA:abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/1:GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA:aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa""".split('\n')

illumina_read1_N = """HWI-6X_9267:1:1:4:1699#ACCACCC/1:TACGGAGGGTGCGAGNGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA:abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB
HWI-6X_9267:1:1:4:390#ACNTCCC/1:GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA:aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa""".split('\n')

illumina_read1_seq_N = """HWI-6X_9267:1:1:4:1699#ACCACCC/1:TACGGAGGGTGCGAGNGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA:abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/1:NACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA:aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa""".split('\n')



illumina_read2 = """HWI-6X_9267:1:1:4:1699#ACCACCC/2:TTTTAAAAAAAAGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTAAAAAAAAACCCCCCCGGGGGGGGTTTTTTTAATTATTC:aaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccBcccccccccccccccc```````BBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/2:ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG:aaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbb""".split('\n')

illumina_read2_N = """HWI-6X_9267:1:1:4:1699#ACCACCC/2:TTTTAAAAAAAAGGGGGNGGGGGCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTAAAAAAAAACCCCCCCGGGGGGGGTTTTTTTAATTATTC:aaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccBcccccccccccccccc```````BBBB
HWI-6X_9267:1:1:4:390#ACNTCCC/2:ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG:aaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbb""".split('\n')

illumina_read2_seq_N = """HWI-6X_9267:1:1:4:1699#ACCACCC/2:TTTTAAAAAAAAGGGGGNGGGGGCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTAAAAAAAAACCCCCCCGGGGGGGGTTTTTTTAATTATTC:aaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccBcccccccccccccccc```````BBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/2:ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG:aaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbb""".split('\n')

barcode_to_sample_id1 = {'GGTGGT':'Samp2', 'GGAGGT':'SAMP_1'}

expected_5prime_seqs_file1 = """>Samp2_0 HWI-6X_9267:1:1:4:1699#ACCACCC
TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
>SAMP_1_1 HWI-6X_9267:1:1:4:390#ACCTCCC
GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAA""".split('\n')

expected_5prime_qual_file1 = """>Samp2_0 HWI-6X_9267:1:1:4:1699#ACCACCC
abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
>SAMP_1_1 HWI-6X_9267:1:1:4:390#ACCTCCC
aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa""".split('\n')

expected_5prime_seqs_file2 = """>Samp2_0 HWI-6X_9267:1:1:4:1699#ACCACCC
TACGGAGGGTGCGAGNGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
>SAMP_1_1 HWI-6X_9267:1:1:4:390#ACCTCCC
NACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAA""".split('\n')

expected_5prime_qual_file2 = """>Samp2_0 HWI-6X_9267:1:1:4:1699#ACCACCC
abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
>SAMP_1_1 HWI-6X_9267:1:1:4:390#ACCTCCC
aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa""".split('\n')

expected_seqs_file1 = """>Samp2_0 HWI-6X_9267:1:1:4:1699#ACCACCC
TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGTTTTTTTTTAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCTTTTTTTTAAAA
>SAMP_1_1 HWI-6X_9267:1:1:4:390#ACCTCCC
GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT""".split('\n')

expected_qual_file1 = """>Samp2_0 HWI-6X_9267:1:1:4:1699#ACCACCC
abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacccccccccccccccccbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
>SAMP_1_1 HWI-6X_9267:1:1:4:390#ACCTCCC
aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaa""".split('\n')

barcode_to_sample_id2 = {'GGTGGT':'Samp2'}

expected_seqs_file2 = """>Samp2_42 HWI-6X_9267:1:1:4:1699#ACCACCC
TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGTTTTTTTTTAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCTTTTTTTTAAAA
>Unassigned_43 HWI-6X_9267:1:1:4:390#ACCTCCC
GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT""".split('\n')

expected_qual_file2 = """>Samp2_42 HWI-6X_9267:1:1:4:1699#ACCACCC
abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacccccccccccccccccbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
>Unassigned_43 HWI-6X_9267:1:1:4:390#ACCTCCC
aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaa""".split('\n')

barcode_to_sample_id3 = barcode_to_sample_id1

expected_seqs_file3 = """>Samp2_0 HWI-6X_9267:1:1:4:1699#ACCACCC
TACGGAGGGTGCGAGNGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGTTTTTTTTTAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGCCCCCNCCCCCTTTTTTTTAAAA
>SAMP_1_1 HWI-6X_9267:1:1:4:390#ACCTCCC
NACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT""".split('\n')

expected_qual_file3 = """>Samp2_0 HWI-6X_9267:1:1:4:1699#ACCACCC
abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacccccccccccccccccbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
>SAMP_1_1 HWI-6X_9267:1:1:4:390#ACCTCCC
aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaa""".split('\n')


barcode_to_sample_id4 = {'ACAGCTAGCTTG':'S1',
                         'CATATCGCAGTT':'S2'}

illumina_read3_bc_in_seq = """HWI-EAS440_0386:1:23:19516:1031#0/1:ACAGCTAGCTTGTACGCAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGGTGTCTTGAGTACAGTAGAGGCAGGCGGAATTCGTGGGG:gggggggeggcffffcGddd\_``_gggggggggggfgggggegggggcgggggggggggeeffafdcdfgbdggdbe]fbfdddddbdadadcddaf`abb`cVNRNUScaa``aOY]]]_[_BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
HWI-EAS440_0386:1:23:19660:1034#0/1:CATATCGCAGTTTACGCAAGGTCCGAGCGTTGTCCGGAATCATTGGGCGTAAAGGGTACGTAGGCGGGTAAGCAAGTTAGAAGTGAAATCCTATAGCTCAACTATAGTAAGCTTTTAAAACTGCTCATCTTGAGGTATGGAAGGGAAAGTGGAATTCCTAGTTA:fhhghhhhfhghhhhcHdcddccddhhhhhhfhhhhghhhdghhhhhhhhhhfhhhdhghgghhhhhhbfdfdbagdgdgfffafa]dad_acdabZcaabad[a__^_`cbddefb_cd^]_L\]U_]^aaZ___]bBBBBBBBBBBBBBBBBBBBBBBBBBB""".split('\n')

expected_seq_file_bc_in_seq1 = """>S1_0 HWI-EAS440_0386:1:23:19516:1031#0
TACGCAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTG
>S2_1 HWI-EAS440_0386:1:23:19660:1034#0
TACGCAAGGTCCGAGCGTTGTCCGGAATCATTGGGCGTAAAGGGTACGTAGGCGGGTAAGCAAGTTAGAAGTGAAATCCTATAGCTCAACTATAGTAAGCTTTTAAAACTGCTCATCTTGAGGTAT""".split('\n')

expected_qual_file_bc_in_seq1 = """>S1_0 HWI-EAS440_0386:1:23:19516:1031#0
fffcGddd\_``_gggggggggggfgggggegggggcgggggggggggeeffafdcdfgbdggdbe]fbfdddddbdadadcddaf`abb`cVNRNUScaa``aOY]]]_[_
>S2_1 HWI-EAS440_0386:1:23:19660:1034#0
hhhcHdcddccddhhhhhhfhhhhghhhdghhhhhhhhhhfhhhdhghgghhhhhhbfdfdbagdgdgfffafa]dad_acdabZcaabad[a__^_`cbddefb_cd^]_L\]U_]^aaZ___]b""".split('\n')

mapping_f = """#SampleID	BarcodeSequence Something Description
Samp2	GGTGGT	True
SAMP_1	GGAGGT	False	Some description
dflsdflsdfsdfsdfsd	GGttAA""".split('\n')

if __name__ == "__main__":
    main()