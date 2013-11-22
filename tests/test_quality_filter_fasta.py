#!/usr/bin/env python
#file test_quality_filter_fasta.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "William Walters"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Development"

from os.path import exists, join
from shutil import rmtree
from StringIO import StringIO

from numpy import array
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from cogent import DNA

from qiime.util import get_tmp_filename, create_dir
from qiime.parse import parse_qual_score
from qiime.quality_filter_fasta import (
    quality_filter_sequences, process_seqs, quality_filter_seq,
    check_primer_mismatch, check_rev_primer_mismatch, check_min_seq_length,
    check_max_seq_length, truncate_ambig_bases, check_ambig_count,
    check_sliding_qual_window, check_average_quality, check_homopolymers,
    check_map, get_ids_primers, format_detailed_data, format_log_data,
    write_qual_line, count_mismatches, expand_degeneracies, make_histograms,
    count_ambig, pair_hmm_align_unaligned_seqs,
    local_align_primer_seq, check_window_qual_scores, seq_exceeds_homopolymers
)

class FakeOutFile(object):
    
    def __init__(self, name="test_file"):
        self.data = ""
        self.name = name
    
    def write(self,s):
        self.data += s
        
class TopLevelTests(TestCase):
    """Tests of top-level functions"""
    
    def setUp(self):
        """ Declares variables, directories for use with unit tests"""
        
        # Input data
        self.valid_fasta_file_no_errors = valid_fasta_file_no_errors
        self.valid_qual_file_no_errors = valid_qual_file_no_errors
        self.valid_mapping_data_golay = valid_mapping_data_golay
        self.valid_fasta_file_no_errors2 = valid_fasta_file_no_errors2
        self.valid_qual_file_no_errors2 = valid_qual_file_no_errors2
        self.fasta_file_ambig_bases = fasta_file_ambig_bases
        self.mapping_data_aligned_primer = mapping_data_aligned_primer
        self.no_matching_sampleid_mapping = no_matching_sampleid_mapping
        self.fasta_file_diff_lengths = fasta_file_diff_lengths
        
        self.valid_fasta_file = get_tmp_filename(prefix = "valid_fasta_file_",
         suffix = ".fasta")
        tmp_file = open(self.valid_fasta_file, 'w')
        tmp_file.write(self.valid_fasta_file_no_errors)
        tmp_file.close()
        
        self.valid_fasta_file2 = get_tmp_filename(prefix = "valid_fasta_file2_",
         suffix = ".fasta")
        tmp_file = open(self.valid_fasta_file2, 'w')
        tmp_file.write(self.valid_fasta_file_no_errors2)
        tmp_file.close()
        
        self.fasta_ambi_bases = get_tmp_filename(prefix = "fasta_ambig_bases_",
         suffix = ".fasta")
        tmp_file = open(self.fasta_ambi_bases, 'w')
        tmp_file.write(self.fasta_file_ambig_bases)
        tmp_file.close()
        
        self.valid_qual_file = get_tmp_filename(prefix = "valid_qual_file_",
         suffix = ".qual")
        tmp_file = open(self.valid_qual_file, 'w')
        tmp_file.write(self.valid_qual_file_no_errors)
        tmp_file.close()
        
        self.valid_qual_file2 = get_tmp_filename(prefix = "valid_qual_file2_",
         suffix = ".qual")
        tmp_file = open(self.valid_qual_file2, 'w')
        tmp_file.write(self.valid_qual_file_no_errors2)
        tmp_file.close()
        
        self.valid_mapping_file = get_tmp_filename(prefix = "valid_mapping_",
         suffix = ".txt")
        tmp_file = open(self.valid_mapping_file, 'w')
        tmp_file.write("\n".join(self.valid_mapping_data_golay))
        tmp_file.close()
        
        self.mapping_file_aligned_primer =\
         get_tmp_filename(prefix = "mapping_aligned_primer_",
         suffix = ".txt")
        tmp_file = open(self.mapping_file_aligned_primer, 'w')
        tmp_file.write("\n".join(self.mapping_data_aligned_primer))
        tmp_file.close()
        
        self.mapping_file_no_sampleids =\
         get_tmp_filename(prefix = "mapping_no_sampleids",
         suffix = ".txt")
        tmp_file = open(self.mapping_file_no_sampleids, 'w')
        tmp_file.write("\n".join(self.no_matching_sampleid_mapping))
        tmp_file.close()
        
        self.fasta_file_diff_lengths_fp =\
         get_tmp_filename(prefix = "fasta_file_diff_lens_",
         suffix = ".txt")
        tmp_file = open(self.fasta_file_diff_lengths_fp, "w")
        tmp_file.write("\n".join(self.fasta_file_diff_lengths))
        tmp_file.close()
        
        self.output_dir = get_tmp_filename(prefix = "quality_filter_fasta_",
         suffix = "/")
        create_dir(self.output_dir)

        self._files_to_remove = \
         [self.valid_fasta_file, self.valid_qual_file,
          self.valid_mapping_file, self.valid_fasta_file2,
          self.valid_qual_file2, self.fasta_ambi_bases,
          self.mapping_file_aligned_primer, self.mapping_file_no_sampleids,
          self.fasta_file_diff_lengths_fp]
          
        # Output data
        self.valid_fasta_demultiplexed = valid_fasta_demultiplexed
        self.log_default_demultiplexed = log_default_demultiplexed
        self.valid_fasta_twofiles_demultiplexed =\
         valid_fasta_twofiles_demultiplexed
        self.log_default_twofile_demultiplexed =\
         log_default_twofile_demultiplexed
        self.log_min_seq_length_no_seqs_written =\
         log_min_seq_length_no_seqs_written
        self.fasta_max_seq_length_changed =\
         fasta_max_seq_length_changed
        self.log_max_seq_length_changed =\
         log_max_seq_length_changed
        self.log_ave_qual_changed = log_ave_qual_changed
        self.fasta_primer_retained = fasta_primer_retained
        self.log_primer_retained = log_primer_retained
        self.output_fasta_ambig_bases = output_fasta_ambig_bases
        self.log_ambig_bases = log_ambig_bases
        self.log_suppress_ambig_check = log_suppress_ambig_check
        self.log_changed_homopolymer = log_changed_homopolymer
        self.log_suppress_homopolymer = log_suppress_homopolymer
        self.fasta_output_one_mm_allowed = fasta_output_one_mm_allowed
        self.log_one_mm_allowed = log_one_mm_allowed
        self.sliding_window_fasta_out = sliding_window_fasta_out
        self.sliding_window_log = sliding_window_log
        self.sliding_window_discard_log = sliding_window_discard_log
        self.fasta_out_rev_primer_trunc = fasta_out_rev_primer_trunc
        self.log_out_rev_primer_trunc = log_out_rev_primer_trunc
        self.fasta_out_rev_primer_remove = fasta_out_rev_primer_remove
        self.log_out_rev_primer_remove = log_out_rev_primer_remove
        self.fasta_out_rev_primer_mm = fasta_out_rev_primer_mm
        self.log_out_rev_primer_mm = log_out_rev_primer_mm
        self.log_out_record_qual_default = log_out_record_qual_default
        self.qual_out_record_default = qual_out_record_default
        self.log_out_trunc_ambi = log_out_trunc_ambi
        self.qual_out_trunc_rev_primer = qual_out_trunc_rev_primer
        self.qual_out_trunc_sliding_window = qual_out_trunc_sliding_window
        self.no_matching_sampleids_fasta = no_matching_sampleids_fasta
        self.no_matching_sampleids_log = no_matching_sampleids_log
        self.sampleids_suppressed_fasta = sampleids_suppressed_fasta
        self.sampleids_suppressed_log = sampleids_suppressed_log
        self.enable_all_tests_main_log = enable_all_tests_main_log
        self.enable_all_tests_detailed_log = enable_all_tests_detailed_log
        self.enable_all_tests_main_log_optional_tests =\
         enable_all_tests_main_log_optional_tests
        self.enable_all_tests_detailed_log_optional_tests =\
         enable_all_tests_detailed_log_optional_tests
        self.expected_fasta_primer_tests_suppressed =\
         expected_fasta_primer_tests_suppressed
        self.expected_log_primer_tests_suppressed =\
         expected_log_primer_tests_suppressed
    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)
            
    def test_quality_filter_sequences(self):
        """ Test overall functionality with single fasta, qual file """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_default_demultiplexed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n', 
         '242\t0\t0\n', '252\t1\t0\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_two_files(self):
        """ Test overall functionality with two fasta, qual files """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_twofiles_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_default_twofile_demultiplexed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t2\t2\n',
         '242\t0\t1\n', '252\t2\t0\n', '262\t2\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_no_qual(self):
        """ Test overall functionality with fasta, no qual files """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], output_dir = self.output_dir)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_default_demultiplexed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n', 
         '242\t0\t0\n', '252\t1\t0\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_seq_lengths(self):
        """ Test overall functionality with altered seq length settings """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         min_seq_len=233, output_dir = self.output_dir)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        # Should be empty, all seqs filtered out
        expected_fasta = []
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_min_seq_length_no_seqs_written
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t1\t0\n',
         '250\t1\t0\n', '260\t1\t0\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_seq_lengths(self):
        """ Handles max seq length changes"""
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, max_seq_len=240)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.fasta_max_seq_length_changed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_max_seq_length_changed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t2\t2\n',
         '242\t0\t0\n', '252\t2\t0\n', '262\t2\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_diff_qual_score(self):
        """ Test overall functionality with different ave qual score """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, min_qual_score=35)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        # Matches output of max sequence length changed 
        expected_fasta = self.fasta_max_seq_length_changed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_ave_qual_changed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t2\t2\n',
         '242\t0\t0\n', '252\t2\t0\n', '262\t2\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_retain_primer(self):
        """ Test overall functionality with primer retained """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, retain_primer=True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.fasta_primer_retained
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_primer_retained
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t1\t1\n',
         '250\t1\t0\n', '260\t1\t0\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_changed_ambig(self):
        """ Test overall functionality with different ambig bases allowed """
        
        # Default setting should have a single fasta seq output
        quality_filter_sequences(self.valid_mapping_file,
         [self.fasta_ambi_bases], [self.valid_qual_file],
         output_dir = self.output_dir)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.output_fasta_ambig_bases
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_default_demultiplexed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n', 
         '242\t0\t0\n', '252\t1\t0\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
        # Lowering setting to 3 max ambiguous bases should write no seqs
        quality_filter_sequences(self.valid_mapping_file,
         [self.fasta_ambi_bases], [self.valid_qual_file],
         output_dir = self.output_dir, max_ambig = 3)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = []
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_ambig_bases
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t1\t0\n',
         '250\t1\t0\n', '260\t1\t0\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_ambig_disabled(self):
        """ Test overall functionality with disabled ambig base check """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.fasta_ambi_bases], [self.valid_qual_file],
         output_dir = self.output_dir, max_ambig = 3, suppress_ambig_check=True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.output_fasta_ambig_bases
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_suppress_ambig_check
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n', 
         '242\t0\t0\n', '252\t1\t0\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_changed_homopolymer(self):
        """ Test overall functionality with different homopolymer len allowed"""
        
        # Should discard all reads with max_homopolymer of 3
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, max_homopolymer = 3)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = []
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_changed_homopolymer
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t1\t0\n',
         '250\t1\t0\n', '260\t1\t0\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_homopolymer_disabled(self):
        """ Test overall functionality with disabled homopolymer check """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, max_homopolymer = 3,
         suppress_homopolymer_check = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_suppress_homopolymer
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n',
         '242\t0\t0\n', '252\t1\t0\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_changed_primer_mismatch(self):
        """ Test overall functionality with different primer mm allowed """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, max_primer_mismatch = 1)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.fasta_output_one_mm_allowed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_one_mm_allowed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n',
         '242\t0\t1\n', '252\t1\t1\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_align_forward_primer(self):
        """ Test overall functionality with aligned forward primer """
        
        quality_filter_sequences(self.mapping_file_aligned_primer,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, local_align_forward_primer = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_default_demultiplexed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n', 
         '242\t0\t0\n', '252\t1\t0\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_suppress_primer_check(self):
        """ Test overall functionality with aligned forward primer """
        
        quality_filter_sequences(self.mapping_file_aligned_primer,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, suppress_primer_check = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        # Should look the same as the input fasta
        expected_fasta = self.expected_fasta_primer_tests_suppressed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.expected_log_primer_tests_suppressed
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t1\t1\n',
         '250\t1\t1\n', '260\t1\t1\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_sliding_qual_window(self):
        """ Test overall functionality with sliding window check enabled"""
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, qual_score_window=10,
         min_qual_score = 32, min_seq_len=50)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.sliding_window_fasta_out
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.sliding_window_log
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '162\t0\t2\n', 
         '172\t0\t0\n', '182\t0\t0\n', '192\t0\t0\n', '202\t0\t0\n',
         '212\t0\t0\n', '222\t0\t0\n', '232\t2\t0\n', '242\t0\t0\n',
         '252\t2\t0\n', '262\t2\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_sliding_window_discard(self):
        """ Test overall functionality with seqs discarded for bad windows """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, qual_score_window=10,
         min_qual_score = 32, min_seq_len=50, discard_bad_windows = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        # Should discard all reads
        actual_fasta = [line for line in output_fasta]
        expected_fasta = []
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.sliding_window_discard_log
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t2\t0\n',
         '250\t2\t0\n', '260\t2\t0\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_reverse_primer_truncate_only(self):
        """ Test overall functionality with reverse primer truncation """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, reverse_primers='truncate_only')
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.fasta_out_rev_primer_trunc
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_out_rev_primer_trunc
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '218\t0\t2\n',
         '228\t0\t0\n', '238\t2\t0\n', '248\t0\t1\n', '258\t4\t0\n',
         '268\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_reverse_primer_remove_option(self):
        """ Test overall functionality with reverse primer remove option """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, reverse_primers='truncate_remove')
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.fasta_out_rev_primer_remove
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_out_rev_primer_remove
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '218\t0\t2\n',
         '228\t0\t0\n', '238\t2\t0\n', '248\t0\t0\n', '258\t4\t0\n',
         '268\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_reverse_primer_mismatches(self):
        """ Test overall functionality with reverse primer mismatch changes """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, reverse_primers='truncate_only',
         reverse_primer_mismatches = 1)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.fasta_out_rev_primer_mm
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_out_rev_primer_mm
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '218\t0\t2\n',
         '228\t0\t1\n', '238\t2\t0\n', '248\t0\t0\n', '258\t4\t0\n',
         '268\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_record_quality_scores(self):
        """ Records quality scores, with and without truncations correctly """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, record_qual_scores = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        output_qual = open(join(self.output_dir, "seqs.qual"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_out_record_qual_default
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n', 
         '242\t0\t0\n', '252\t1\t0\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        actual_qual = [line for line in output_qual]
        expected_qual = self.qual_out_record_default
        self.assertEqual(actual_qual, expected_qual)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        output_qual.close()
        
    def test_quality_filter_sequences_truncate_ambig_bases(self):
        """ Test overall functionality with ambiguous base truncation """
        
        # Default setting should have a single fasta seq output
        # Also test qual scores recording with truncation
        quality_filter_sequences(self.valid_mapping_file,
         [self.fasta_ambi_bases], [self.valid_qual_file],
         output_dir = self.output_dir, truncate_ambi_bases = True,
         min_seq_len = 10, record_qual_scores = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        output_qual = open(join(self.output_dir, "seqs.qual"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n',
         'CTCCCGTAGGAGTC\n']
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.log_out_trunc_ambi
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '14\t0\t1\n', '24\t0\t0\n',
         '34\t0\t0\n', '44\t0\t0\n', '54\t0\t0\n', '64\t0\t0\n', '74\t0\t0\n',
         '84\t0\t0\n', '94\t0\t0\n', '104\t0\t0\n', '114\t0\t0\n',
         '124\t0\t0\n', '134\t0\t0\n', '144\t0\t0\n', '154\t0\t0\n',
         '164\t0\t0\n', '174\t0\t0\n', '184\t0\t0\n', '194\t0\t0\n',
         '204\t0\t0\n', '214\t0\t0\n', '224\t0\t0\n', '234\t1\t0\n',
         '244\t0\t0\n', '254\t1\t0\n', '264\t1\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        actual_qual = [line for line in output_qual]
        expected_qual = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n',
         '40 40 40 40 40 40 40 40 38 38 39 40 40 40\n']
        self.assertEqual(actual_qual, expected_qual)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        output_qual.close()
        
        # Should have no output if min length left at default, too short
        # after truncation
        quality_filter_sequences(self.valid_mapping_file,
         [self.fasta_ambi_bases], [self.valid_qual_file],
         output_dir = self.output_dir, truncate_ambi_bases = True,
         record_qual_scores = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        output_qual = open(join(self.output_dir, "seqs.qual"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = []
        self.assertEqual(actual_fasta, expected_fasta)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t1\t0\n',
         '250\t1\t0\n', '260\t1\t0\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        actual_qual = [line for line in output_qual]
        expected_qual = []
        self.assertEqual(actual_qual, expected_qual)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        output_qual.close()
        
    def test_quality_filter_sequences_reverse_primer_record_qual(self):
        """ Test qual score recording with reverse primer truncation """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, reverse_primers='truncate_only',
         record_qual_scores = True)
        
        output_qual = open(join(self.output_dir, "seqs.qual"), "U")
        
        actual_qual = [line for line in output_qual]
        expected_qual = self.qual_out_trunc_rev_primer
        self.assertEqual(actual_qual, expected_qual)
        output_qual.close()
        
    def test_quality_filter_sequences_sliding_qual_window_record_qual(self):
        """ Records qual scores with sliding window check enabled"""
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file, self.valid_fasta_file2],
         [self.valid_qual_file, self.valid_qual_file2],
         output_dir = self.output_dir, qual_score_window=10,
         min_qual_score = 32, min_seq_len=50, record_qual_scores = True)
        
        output_qual = open(join(self.output_dir, "seqs.qual"), "U")
        
        actual_qual = [line for line in output_qual]
        expected_qual = self.qual_out_trunc_sliding_window
        self.assertEqual(actual_qual, expected_qual)
        output_qual.close()

        
    def test_quality_filter_sequences_sampleid_suppression(self):
        """ Test overall functionality with sampleid checks suppressed """
        
        # Without matching SampleIDs, should fail to detect primers
        quality_filter_sequences(self.mapping_file_no_sampleids,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.no_matching_sampleids_fasta
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.no_matching_sampleids_log
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t1\t1\n',
         '250\t1\t1\n', '260\t1\t1\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
        # Should remove primers with suppress sampleID enabled
        quality_filter_sequences(self.mapping_file_no_sampleids,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, suppress_sampleid_check = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        
        # Primers removed from all reads, all written.
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.sampleids_suppressed_fasta
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.sampleids_suppressed_log
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n',
         '242\t0\t1\n', '252\t1\t1\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        
    def test_quality_filter_sequences_enable_all_tests(self):
        """ Test overall functionality with all tests enabled """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, record_qual_scores = True,
         enable_all_checks = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        output_qual = open(join(self.output_dir, "seqs.qual"), "U")
        output_detailed = open(join(self.output_dir,
         "detailed_quality_report.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.enable_all_tests_main_log
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '232\t1\t1\n', 
         '242\t0\t0\n', '252\t1\t0\n', '262\t1\t0\n', '272\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        actual_qual = [line for line in output_qual]
        expected_qual = self.qual_out_record_default
        self.assertEqual(actual_qual, expected_qual)
        
        actual_detailed_log = [line for line in output_detailed]
        expected_detailed_log = self.enable_all_tests_detailed_log
        self.assertEqual(actual_detailed_log, expected_detailed_log)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        output_qual.close()
        output_detailed.close()
        
    def test_quality_filter_sequences_enable_all_optional_tests(self):
        """ Test overall functionality with optional tests enabled """
        
        quality_filter_sequences(self.valid_mapping_file,
         [self.valid_fasta_file], [self.valid_qual_file],
         output_dir = self.output_dir, record_qual_scores = True,
         enable_all_checks = True, min_qual_score = 35, max_ambig=3,
         qual_score_window = 10, reverse_primers = 'truncate_only',
         truncate_ambi_bases = True)
        
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        output_log = open(join(self.output_dir, "quality_filter_log.txt"), "U")
        output_hist = open(join(self.output_dir, "histograms.txt"), "U")
        output_qual = open(join(self.output_dir, "seqs.qual"), "U")
        output_detailed = open(join(self.output_dir,
         "detailed_quality_report.txt"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = []
        self.assertEqual(actual_fasta, expected_fasta)
        
        # Slice off the first 5 lines, as contains dynamic data regarding
        # temp filepaths used
        actual_log = [line for line in output_log][5:]
        expected_log = self.enable_all_tests_main_log_optional_tests
        self.assertEqual(actual_log, expected_log)
        
        actual_histogram = [line for line in output_hist]
        expected_hist = ['Length\tBefore\tAfter\n', '240\t1\t0\n',
         '250\t1\t0\n', '260\t1\t0\n', '270\t0\t0']
        self.assertEqual(actual_histogram, expected_hist)
        
        actual_qual = [line for line in output_qual]
        expected_qual = []
        self.assertEqual(actual_qual, expected_qual)
        
        actual_detailed_log = [line for line in output_detailed]
        expected_detailed_log =\
         self.enable_all_tests_detailed_log_optional_tests
        self.assertEqual(actual_detailed_log, expected_detailed_log)
        
        output_fasta.close()
        output_log.close()
        output_hist.close()
        output_qual.close()
        output_detailed.close()
        
    def test_process_seqs_with_qual_data(self):
        """ process_seqs handles fasta/qual file """
        
        ids_primers = {'PC.634':['CATGCTTC'], 'PC.481':['CATGCTGC', 'CATGTTGC']}
        ids_rev_primers = {}
        process_seqs([self.valid_fasta_file], [self.valid_qual_file],
         ids_primers, ids_rev_primers, output_dir = self.output_dir)
         
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        output_fasta.close()
        
    def test_process_seqs_with_no_qual_data(self):
        """ process_seqs handles fasta file alone """
        
        ids_primers = {'PC.634':['CATGCTTC'], 'PC.481':['CATGCTGC', 'CATGTTGC']}
        ids_rev_primers = {}
        process_seqs([self.valid_fasta_file], [],
         ids_primers, ids_rev_primers, output_dir = self.output_dir)
         
        output_fasta = open(join(self.output_dir, "seqs.fna"), "U")
        
        actual_fasta = [line for line in output_fasta]
        expected_fasta = self.valid_fasta_demultiplexed
        self.assertEqual(actual_fasta, expected_fasta)
        
        output_fasta.close()
        
    def test_process_seqs_errors_on_mismatched_labels(self):
        """ process_seqs errors when fasta/qual labels are not matched """
        
        ids_primers = {'PC.634':['CATGCTTC'], 'PC.481':['CATGCTGC', 'CATGTTGC']}
        ids_rev_primers = {}
        self.assertRaises(ValueError, process_seqs, [self.valid_fasta_file],
         [self.valid_qual_file2], ids_primers, ids_rev_primers,
         output_dir = self.output_dir)

        
    def test_process_seqs_errors_on_diff_seq_lens(self):
        """ process_seqs errors when fasta/qual seq lens are not matched """
        
        ids_primers = {'PC.634':['CATGCTTC'], 'PC.481':['CATGCTGC', 'CATGTTGC']}
        ids_rev_primers = {}
        self.assertRaises(ValueError, process_seqs,
         [self.fasta_file_diff_lengths_fp], [self.valid_qual_file],
         ids_primers, ids_rev_primers, output_dir = self.output_dir)

    def test_quality_filter_seq(self):
        """ quality_filter_seq filters seq, logs data """
        
        ids_primers = {'PC.634':['CATGCTTC'], 'PC.481':['CATGCTGC', 'CATGTTGC']}
        ids_rev_primers = {}
        final_log_data = {'seq_counts':0,
                      'seqs_written':0,
                      'fasta_files':'fasta_files',
                      'qual_files':'qual_files',
                      'below_min_seq_len':0,
                      'max_seq_len':0,
                      'min_ave_qual_score':0,
                      'max_ambig':0,
                      'ambig_base_truncation':0,
                      'too_short_after_ambig_truncation':0,
                      'max_homopolymer':0,
                      'max_primer_mismatch':0,
                      'low_qual_window_found':0,
                      'too_short_after_window_truncation':0,
                      'low_qual_window_discarded':0,
                      'rev_primers_found':0,
                      'max_rev_primer_mismatch':0,
                      'seqs_discarded_no_rev_primer':0,
                      'too_short_after_revprimer_truncation':0,
                      'seq_ids_not_in_mapping':[]
                      }
        detailed_quality_data = {'seq_order':[]}
        detailed_quality_data['all_tests'] = ['below_min_seq_len',
         'exceeds_max_seq_len', 'exceeds_max_primer_mismatch',
         'rev_primer_found', 'exceeds_max_rev_primer_mismatch',
         'seqs_discarded_for_no_rev_primer', 'too_short_after_revprimer_truncation',
         'truncation_for_ambig_base', 'too_short_after_ambig_truncation',
         'exceeds_max_ambig', 'exceeds_max_homopolymer', 'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping']
        fasta_seq, qual_seq, final_log_data, detailed_quality_data =\
         quality_filter_seq(fasta_label = ">test1", fasta_seq="ACTACGGATACG",
         qual_seq=array([30, 30, 25, 25, 15, 25, 30, 35, 25, 27, 30, 31]),
         final_log_data = final_log_data,
         detailed_quality_data = detailed_quality_data,
         ids_primers = ids_primers, ids_rev_primers = ids_rev_primers,
         min_seq_len = 10)
         
        expected_fasta_seq = "ACTACGGATACG"
        expected_qual_seq = array([30, 30, 25, 25, 15, 25, 30, 35, 25,
         27, 30, 31])
        expected_final_log_data = {'rev_primers_found': 0,
         'max_homopolymer': 0, 'low_qual_window_found': 0, 'seq_counts': 0,
         'max_seq_len': 0, 'low_qual_window_discarded': 0,
         'ambig_base_truncation': 0, 'seqs_discarded_no_rev_primer': 0,
         'max_ambig': 0, 'qual_files': 'qual_files',
         'fasta_files': 'fasta_files', 'seq_ids_not_in_mapping': ['>test1'],
         'max_primer_mismatch': 0, 'max_rev_primer_mismatch': 0,
         'too_short_after_revprimer_truncation': 0, 'seqs_written': 0,
         'too_short_after_window_truncation': 0,
         'too_short_after_ambig_truncation': 0, 'min_ave_qual_score': 0,
         'below_min_seq_len': 0}
        expected_detailed_quality_data = {'seq_order': [],
         'all_tests': ['below_min_seq_len', 'exceeds_max_seq_len',
         'exceeds_max_primer_mismatch', 'rev_primer_found',
         'exceeds_max_rev_primer_mismatch', 'seqs_discarded_for_no_rev_primer',
         'too_short_after_revprimer_truncation', 'truncation_for_ambig_base',
         'too_short_after_ambig_truncation', 'exceeds_max_ambig',
         'exceeds_max_homopolymer', 'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping'],
         '>test1': {}}
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_quality_filter_seq_all_tests(self):
        """ quality_filter_seq does all tests when enabled """
        
        ids_primers = {'all_primers':['CATGCTTC', 'CATGCTGC', 'CATGTTGC']}
        ids_rev_primers = {}
        final_log_data = {'seq_counts':0,
                      'seqs_written':0,
                      'fasta_files':'fasta_files',
                      'qual_files':'qual_files',
                      'below_min_seq_len':0,
                      'max_seq_len':0,
                      'min_ave_qual_score':0,
                      'max_ambig':0,
                      'ambig_base_truncation':0,
                      'too_short_after_ambig_truncation':0,
                      'max_homopolymer':0,
                      'max_primer_mismatch':0,
                      'low_qual_window_found':0,
                      'too_short_after_window_truncation':0,
                      'low_qual_window_discarded':0,
                      'rev_primers_found':0,
                      'max_rev_primer_mismatch':0,
                      'seqs_discarded_no_rev_primer':0,
                      'too_short_after_revprimer_truncation':0,
                      'seq_ids_not_in_mapping':[]
                      }
        detailed_quality_data = {'seq_order':[]}
        detailed_quality_data['all_tests'] = ['below_min_seq_len',
         'exceeds_max_seq_len', 'exceeds_max_primer_mismatch',
         'rev_primer_found', 'exceeds_max_rev_primer_mismatch',
         'seqs_discarded_for_no_rev_primer', 'too_short_after_revprimer_truncation',
         'truncation_for_ambig_base', 'too_short_after_ambig_truncation',
         'exceeds_max_ambig', 'exceeds_max_homopolymer', 'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping']
        fasta_seq, qual_seq, final_log_data, detailed_quality_data =\
         quality_filter_seq(fasta_label = ">test1", fasta_seq="ACTACGGATACG",
         qual_seq=array([30, 30, 25, 25, 15, 25, 30, 35, 25, 27, 30, 31]),
         final_log_data = final_log_data,
         detailed_quality_data = detailed_quality_data,
         ids_primers = ids_primers, ids_rev_primers = ids_rev_primers,
         min_seq_len = 14, suppress_sampleid_check = True,
         enable_all_checks = True)
         
        expected_fasta_seq = ""
        expected_qual_seq = ""
        expected_final_log_data = {'rev_primers_found': 0,
         'max_homopolymer': 0, 'low_qual_window_found': 0, 'seq_counts': 0,
         'max_seq_len': 0, 'low_qual_window_discarded': 0,
         'ambig_base_truncation': 0, 'seqs_discarded_no_rev_primer': 0,
         'max_ambig': 0, 'qual_files': 'qual_files',
         'fasta_files': 'fasta_files', 'seq_ids_not_in_mapping': [],
         'max_primer_mismatch': 1, 'max_rev_primer_mismatch': 0,
         'too_short_after_revprimer_truncation': 0, 'seqs_written': 0,
         'too_short_after_window_truncation': 0,
         'too_short_after_ambig_truncation': 0, 'min_ave_qual_score': 0,
         'below_min_seq_len': 1}
        expected_detailed_quality_data = {'seq_order': ['>test1'],
         'all_tests': ['below_min_seq_len', 'exceeds_max_seq_len',
         'exceeds_max_primer_mismatch', 'rev_primer_found',
         'exceeds_max_rev_primer_mismatch', 'seqs_discarded_for_no_rev_primer',
         'too_short_after_revprimer_truncation', 'truncation_for_ambig_base',
         'too_short_after_ambig_truncation', 'exceeds_max_ambig',
         'exceeds_max_homopolymer', 'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping'],
         '>test1': {'exceeds_max_homopolymer': 0,
         'below_min_ave_qual_score': 0, 'exceeds_max_ambig': 0,
         'exceeds_max_primer_mismatch': 1, 'exceeds_max_seq_len': 0,
         'below_min_seq_len': 1}}
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_check_primer_mismatch(self):
        """ check_primer_mismatch detects forward primer mismatches """
        
        fasta_label = "test0_1"
        fasta_seq = "AACCTTGG"
        qual_seq = [1,2,3,4,5,6,7,8]
        ids_primers = {'test0':['AACCTTGG'], 'PC.481':['CATGCTGC', 'CATGTTGC']}
        final_log_data = {'seq_counts':0,
                      'seqs_written':0,
                      'fasta_files':'fasta_files',
                      'qual_files':'qual_files',
                      'below_min_seq_len':0,
                      'max_seq_len':0,
                      'min_ave_qual_score':0,
                      'max_ambig':0,
                      'ambig_base_truncation':0,
                      'too_short_after_ambig_truncation':0,
                      'max_homopolymer':0,
                      'max_primer_mismatch':0,
                      'low_qual_window_found':0,
                      'too_short_after_window_truncation':0,
                      'low_qual_window_discarded':0,
                      'rev_primers_found':0,
                      'max_rev_primer_mismatch':0,
                      'seqs_discarded_no_rev_primer':0,
                      'too_short_after_revprimer_truncation':0,
                      'seq_ids_not_in_mapping':[]
                      }
        detailed_quality_data = {'seq_order':[]}
        detailed_quality_data['all_tests'] = ['below_min_seq_len',
         'exceeds_max_seq_len', 'exceeds_max_primer_mismatch',
         'rev_primer_found', 'exceeds_max_rev_primer_mismatch',
         'seqs_discarded_for_no_rev_primer',
         'too_short_after_revprimer_truncation',
         'truncation_for_ambig_base', 'too_short_after_ambig_truncation',
         'exceeds_max_ambig', 'exceeds_max_homopolymer',
         'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping']
        
        failed, fasta_seq, qual_seq, final_log_data =\
         check_primer_mismatch(fasta_label, fasta_seq, qual_seq, final_log_data,
         detailed_quality_data, ids_primers, retain_primer = False,
         max_primer_mismatch = 0, suppress_sampleid_check = False,
         enable_all_checks = False, local_align_forward_primer = False)
         
        expected_failed = False
        expected_fasta_seq = ""
        expected_qual_seq = []
        expected_final_log_data = {'rev_primers_found': 0, 'max_homopolymer': 0,
         'low_qual_window_found': 0, 'seq_counts': 0, 'max_seq_len': 0,
         'low_qual_window_discarded': 0, 'ambig_base_truncation': 0,
         'seqs_discarded_no_rev_primer': 0, 'max_ambig': 0,
         'qual_files': 'qual_files', 'fasta_files': 'fasta_files',
         'seq_ids_not_in_mapping': [], 'max_primer_mismatch': 0,
         'max_rev_primer_mismatch': 0,
         'too_short_after_revprimer_truncation': 0, 'seqs_written': 0,
         'too_short_after_window_truncation': 0,
         'too_short_after_ambig_truncation': 0, 'min_ave_qual_score': 0,
         'below_min_seq_len': 0}

        self.assertEqual(failed, expected_failed)
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        
    def test_check_primer_mismatch_no_sample_id_check(self):
        """ check_primer_mismatch handles case of suppressed sampleid checks """
        
        fasta_label = ">test0_1"
        fasta_seq = "AACCTTGG"
        qual_seq = [1,2,3,4,5,6,7,8]
        ids_primers = {'all_primers':['AACCTTGG', 'CATGCTGC', 'CATGTTGC']}
        final_log_data = {'seq_counts':0,
                      'seqs_written':0,
                      'fasta_files':'fasta_files',
                      'qual_files':'qual_files',
                      'below_min_seq_len':0,
                      'max_seq_len':0,
                      'min_ave_qual_score':0,
                      'max_ambig':0,
                      'ambig_base_truncation':0,
                      'too_short_after_ambig_truncation':0,
                      'max_homopolymer':0,
                      'max_primer_mismatch':0,
                      'low_qual_window_found':0,
                      'too_short_after_window_truncation':0,
                      'low_qual_window_discarded':0,
                      'rev_primers_found':0,
                      'max_rev_primer_mismatch':0,
                      'seqs_discarded_no_rev_primer':0,
                      'too_short_after_revprimer_truncation':0,
                      'seq_ids_not_in_mapping':[]
                      }
        detailed_quality_data = {'seq_order':[]}
        detailed_quality_data['all_tests'] = ['below_min_seq_len',
         'exceeds_max_seq_len', 'exceeds_max_primer_mismatch',
         'rev_primer_found', 'exceeds_max_rev_primer_mismatch',
         'seqs_discarded_for_no_rev_primer',
         'too_short_after_revprimer_truncation',
         'truncation_for_ambig_base', 'too_short_after_ambig_truncation',
         'exceeds_max_ambig', 'exceeds_max_homopolymer',
         'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping']
        
        failed, fasta_seq, qual_seq, final_log_data =\
         check_primer_mismatch(fasta_label, fasta_seq, qual_seq, final_log_data,
         detailed_quality_data, ids_primers, retain_primer = False,
         max_primer_mismatch = 0, suppress_sampleid_check = True,
         enable_all_checks = False, local_align_forward_primer = False)
         
        expected_failed = False
        expected_fasta_seq = ""
        expected_qual_seq = []
        expected_final_log_data = {'rev_primers_found': 0, 'max_homopolymer': 0,
         'low_qual_window_found': 0, 'seq_counts': 0, 'max_seq_len': 0,
         'low_qual_window_discarded': 0, 'ambig_base_truncation': 0,
         'seqs_discarded_no_rev_primer': 0, 'max_ambig': 0,
         'qual_files': 'qual_files', 'fasta_files': 'fasta_files',
         'seq_ids_not_in_mapping': [], 'max_primer_mismatch': 0,
         'max_rev_primer_mismatch': 0,
         'too_short_after_revprimer_truncation': 0, 'seqs_written': 0,
         'too_short_after_window_truncation': 0,
         'too_short_after_ambig_truncation': 0, 'min_ave_qual_score': 0,
         'below_min_seq_len': 0}

        self.assertEqual(failed, expected_failed)
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        
    def test_check_rev_primer_mismatch(self):
        """ check_rev_primer_mismatch detects rev primer index and mismatches"""
        
        fasta_label = "test0_1"
        fasta_seq = "AACCTTGG"
        qual_seq = [1,2,3,4,5,6,7,8]
        ids_rev_primers = {'test0':['CTTG'], 'PC.481':['CATGCTGC', 'CATGTTGC']}
        
        final_log_data = {'seq_counts':0,
                      'seqs_written':0,
                      'fasta_files':'fasta_files',
                      'qual_files':'qual_files',
                      'below_min_seq_len':0,
                      'max_seq_len':0,
                      'min_ave_qual_score':0,
                      'max_ambig':0,
                      'ambig_base_truncation':0,
                      'too_short_after_ambig_truncation':0,
                      'max_homopolymer':0,
                      'max_primer_mismatch':0,
                      'low_qual_window_found':0,
                      'too_short_after_window_truncation':0,
                      'low_qual_window_discarded':0,
                      'rev_primers_found':0,
                      'max_rev_primer_mismatch':0,
                      'seqs_discarded_no_rev_primer':0,
                      'too_short_after_revprimer_truncation':0,
                      'seq_ids_not_in_mapping':[]
                      }
        detailed_quality_data = {'seq_order':[]}
        detailed_quality_data['all_tests'] = ['below_min_seq_len',
         'exceeds_max_seq_len', 'exceeds_max_primer_mismatch',
         'rev_primer_found', 'exceeds_max_rev_primer_mismatch',
         'seqs_discarded_for_no_rev_primer',
         'too_short_after_revprimer_truncation',
         'truncation_for_ambig_base', 'too_short_after_ambig_truncation',
         'exceeds_max_ambig', 'exceeds_max_homopolymer',
         'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping']
        
        failed, fasta_seq, qual_seq, final_log_data =\
         check_rev_primer_mismatch(fasta_label,
                              fasta_seq,
                              qual_seq,
                              final_log_data,
                              detailed_quality_data,
                              ids_rev_primers,
                              reverse_primer_mismatches = 0,
                              reverse_primers = "truncate_only",
                              suppress_sampleid_check = False,
                              enable_all_checks = False)
                              
        expected_failed = False
        expected_fasta_seq = "AAC"
        expected_qual_seq = [1,2,3]
        expected_final_log_data = {'rev_primers_found': 1, 'max_homopolymer': 0,
         'low_qual_window_found': 0, 'seq_counts': 0, 'max_seq_len': 0,
         'low_qual_window_discarded': 0, 'ambig_base_truncation': 0,
         'seqs_discarded_no_rev_primer': 0, 'max_ambig': 0,
         'qual_files': 'qual_files', 'fasta_files': 'fasta_files',
         'seq_ids_not_in_mapping': [], 'max_primer_mismatch': 0,
         'max_rev_primer_mismatch': 0,
         'too_short_after_revprimer_truncation': 0, 'seqs_written': 0,
         'too_short_after_window_truncation': 0,
         'too_short_after_ambig_truncation': 0, 'min_ave_qual_score': 0,
         'below_min_seq_len': 0}

        self.assertEqual(failed, expected_failed)
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        
    def test_check_rev_primer_mismatch_no_sample_id_check(self):
        """ check_rev_primer_mismatch handles suppressed sample id check """
        
        fasta_label = "test0_1"
        fasta_seq = "AACCTTGG"
        qual_seq = [1,2,3,4,5,6,7,8]
        ids_rev_primers = {'rev_primers':['CTTG', 'CATGCTGC', 'CATGTTGC']}
        
        final_log_data = {'seq_counts':0,
                      'seqs_written':0,
                      'fasta_files':'fasta_files',
                      'qual_files':'qual_files',
                      'below_min_seq_len':0,
                      'max_seq_len':0,
                      'min_ave_qual_score':0,
                      'max_ambig':0,
                      'ambig_base_truncation':0,
                      'too_short_after_ambig_truncation':0,
                      'max_homopolymer':0,
                      'max_primer_mismatch':0,
                      'low_qual_window_found':0,
                      'too_short_after_window_truncation':0,
                      'low_qual_window_discarded':0,
                      'rev_primers_found':0,
                      'max_rev_primer_mismatch':0,
                      'seqs_discarded_no_rev_primer':0,
                      'too_short_after_revprimer_truncation':0,
                      'seq_ids_not_in_mapping':[]
                      }
        detailed_quality_data = {'seq_order':[]}
        detailed_quality_data['all_tests'] = ['below_min_seq_len',
         'exceeds_max_seq_len', 'exceeds_max_primer_mismatch',
         'rev_primer_found', 'exceeds_max_rev_primer_mismatch',
         'seqs_discarded_for_no_rev_primer',
         'too_short_after_revprimer_truncation',
         'truncation_for_ambig_base', 'too_short_after_ambig_truncation',
         'exceeds_max_ambig', 'exceeds_max_homopolymer',
         'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping']
        
        failed, fasta_seq, qual_seq, final_log_data =\
         check_rev_primer_mismatch(fasta_label,
                              fasta_seq,
                              qual_seq,
                              final_log_data,
                              detailed_quality_data,
                              ids_rev_primers,
                              reverse_primer_mismatches = 0,
                              reverse_primers = "truncate_only",
                              suppress_sampleid_check = True,
                              enable_all_checks = False)
                              
        expected_failed = False
        expected_fasta_seq = "AAC"
        expected_qual_seq = [1,2,3]
        expected_final_log_data = {'rev_primers_found': 1, 'max_homopolymer': 0,
         'low_qual_window_found': 0, 'seq_counts': 0, 'max_seq_len': 0,
         'low_qual_window_discarded': 0, 'ambig_base_truncation': 0,
         'seqs_discarded_no_rev_primer': 0, 'max_ambig': 0,
         'qual_files': 'qual_files', 'fasta_files': 'fasta_files',
         'seq_ids_not_in_mapping': [], 'max_primer_mismatch': 0,
         'max_rev_primer_mismatch': 0,
         'too_short_after_revprimer_truncation': 0, 'seqs_written': 0,
         'too_short_after_window_truncation': 0,
         'too_short_after_ambig_truncation': 0, 'min_ave_qual_score': 0,
         'below_min_seq_len': 0}

        self.assertEqual(failed, expected_failed)
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        
    def test_check_rev_primer_mismatch_truncate_remove(self):
        """ check_rev_primer_mismatch removes seqs with truncate_remove """
        
        fasta_label = "test0_1"
        fasta_seq = "AACCTTGG"
        qual_seq = [1,2,3,4,5,6,7,8]
        ids_rev_primers = {'test0':['CATGCTGC'],
         'PC.481':['CATGCTGC', 'CATGTTGC']}
        
        final_log_data = {'seq_counts':0,
                      'seqs_written':0,
                      'fasta_files':'fasta_files',
                      'qual_files':'qual_files',
                      'below_min_seq_len':0,
                      'max_seq_len':0,
                      'min_ave_qual_score':0,
                      'max_ambig':0,
                      'ambig_base_truncation':0,
                      'too_short_after_ambig_truncation':0,
                      'max_homopolymer':0,
                      'max_primer_mismatch':0,
                      'low_qual_window_found':0,
                      'too_short_after_window_truncation':0,
                      'low_qual_window_discarded':0,
                      'rev_primers_found':0,
                      'max_rev_primer_mismatch':0,
                      'seqs_discarded_no_rev_primer':0,
                      'too_short_after_revprimer_truncation':0,
                      'seq_ids_not_in_mapping':[]
                      }
        detailed_quality_data = {'seq_order':[]}
        detailed_quality_data['all_tests'] = ['below_min_seq_len',
         'exceeds_max_seq_len', 'exceeds_max_primer_mismatch',
         'rev_primer_found', 'exceeds_max_rev_primer_mismatch',
         'seqs_discarded_for_no_rev_primer',
         'too_short_after_revprimer_truncation',
         'truncation_for_ambig_base', 'too_short_after_ambig_truncation',
         'exceeds_max_ambig', 'exceeds_max_homopolymer',
         'below_min_ave_qual_score',
         'low_qual_window_found', 'discarded_for_low_qual_window',
         'too_short_after_window_truncation', 'seq_id_not_in_mapping']
        
        failed, fasta_seq, qual_seq, final_log_data =\
         check_rev_primer_mismatch(fasta_label,
                              fasta_seq,
                              qual_seq,
                              final_log_data,
                              detailed_quality_data,
                              ids_rev_primers,
                              reverse_primer_mismatches = 0,
                              reverse_primers = "truncate_remove",
                              suppress_sampleid_check = False,
                              enable_all_checks = False)
                              
        expected_failed = True
        expected_fasta_seq = "AACCTTGG"
        expected_qual_seq = [1, 2, 3, 4, 5, 6, 7, 8]
        expected_final_log_data = {'rev_primers_found': 0, 'max_homopolymer': 0,
         'low_qual_window_found': 0, 'seq_counts': 0, 'max_seq_len': 0,
         'low_qual_window_discarded': 0, 'ambig_base_truncation': 0,
         'seqs_discarded_no_rev_primer': 1, 'max_ambig': 0,
         'qual_files': 'qual_files', 'fasta_files': 'fasta_files',
         'seq_ids_not_in_mapping': [], 'max_primer_mismatch': 0,
         'max_rev_primer_mismatch': 1,
         'too_short_after_revprimer_truncation': 0, 'seqs_written': 0,
         'too_short_after_window_truncation': 0,
         'too_short_after_ambig_truncation': 0, 'min_ave_qual_score': 0,
         'below_min_seq_len': 0}

        self.assertEqual(failed, expected_failed)
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        
    def test_check_min_seq_length(self):
        """ check_min_seq_length properly flags seqs below minimum length """
        
        fasta_seq = "AACCTTGG"
        final_log_data = {'below_min_seq_len':10}
        detailed_quality_data = {}
        min_seq_len = 8
        failed, final_log_data = check_min_seq_length(fasta_seq,
                         final_log_data,
                         detailed_quality_data,
                         min_seq_len,
                         enable_all_checks = False,
                         curr_log_key = 'below_min_seq_len')
        expected_failed = False
        expected_final_log_data = {'below_min_seq_len': 10}
        expected_detailed_quality_data = {}
        
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_check_min_seq_length_fails(self):
        """ check_min_seq_length properly flags seqs below minimum length """
        
        fasta_seq = "AACCTTGG"
        final_log_data = {'below_min_seq_len':10}
        detailed_quality_data = {}
        min_seq_len = 10
        failed, final_log_data = check_min_seq_length(fasta_seq,
                         final_log_data,
                         detailed_quality_data,
                         min_seq_len,
                         enable_all_checks = True,
                         curr_log_key = 'below_min_seq_len')
        expected_failed = True
        expected_final_log_data = {'below_min_seq_len': 11}
        expected_detailed_quality_data = {'below_min_seq_len': 1}
        
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data) 
        
    def test_check_max_seq_length(self):
        """ check_max_seq_length properly flags seqs above maximum length """
        
        fasta_seq = "AACCTTGG"
        final_log_data = {'max_seq_len':10}
        detailed_quality_data = {}
        max_seq_len = 8
        failed, final_log_data = check_max_seq_length(fasta_seq,
                         final_log_data,
                         detailed_quality_data,
                         max_seq_len,
                         enable_all_checks = False)
        expected_failed = False
        expected_final_log_data = {'max_seq_len': 10}
        expected_detailed_quality_data = {}
        
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_check_max_seq_length_fails(self):
        """ check_max_seq_length properly flags seqs above maximum length """
        
        fasta_seq = "AACCTTGG"
        final_log_data = {'max_seq_len':10}
        detailed_quality_data = {}
        max_seq_len = 6
        failed, final_log_data = check_max_seq_length(fasta_seq,
                         final_log_data,
                         detailed_quality_data,
                         max_seq_len,
                         enable_all_checks = True)
        expected_failed = True
        expected_final_log_data = {'max_seq_len': 11}
        expected_detailed_quality_data = {'exceeds_max_seq_len':1}
        
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_truncate_ambig_bases(self):
        """ truncate_ambig_bases properly truncates at first ambiguous base """
        
        fasta_seq = "ACTCGATASAT"
        qual_seq = [1,2,3,4,5,6,7,8,9,10,11]
        final_log_data = {'ambig_base_truncation':0}
        detailed_quality_data = {}
        
        fasta_seq, qual_seq, final_log_data, truncated =\
         truncate_ambig_bases(fasta_seq,
                         qual_seq,
                         final_log_data,
                         detailed_quality_data,
                         enable_all_checks = False)
        
        expected_fasta_seq = "ACTCGATA"
        expected_qual_seq = [1,2,3,4,5,6,7,8]
        expected_final_log_data = {'ambig_base_truncation':1}
        expected_detailed_quality_data = {}
        
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_truncate_ambig_bases_detailed(self):
        """ truncate_ambig_bases properly truncates at first ambiguous base """
        
        fasta_seq = "ACTCGATASAT"
        qual_seq = [1,2,3,4,5,6,7,8,9,10,11]
        final_log_data = {'ambig_base_truncation':0}
        detailed_quality_data = {}
        
        fasta_seq, qual_seq, final_log_data, truncated =\
         truncate_ambig_bases(fasta_seq,
                         qual_seq,
                         final_log_data,
                         detailed_quality_data,
                         enable_all_checks = True)
        
        expected_fasta_seq = "ACTCGATA"
        expected_qual_seq = [1,2,3,4,5,6,7,8]
        expected_final_log_data = {'ambig_base_truncation':1}
        expected_detailed_quality_data = {'truncation_for_ambig_base':1}
        
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_check_ambig_count(self):
        """ check_ambig_count properly flags for ambiguous base count """
        
        fasta_seq = "ACTCGATASAT"
        final_log_data = {'max_ambig':0}
        detailed_quality_data = {}
        
        failed, final_log_data =\
         check_ambig_count(fasta_seq,
                           final_log_data,
                           detailed_quality_data,
                           max_ambig = 6,
                           enable_all_checks = False)
                           
        expected_failed = False
        expected_final_log_data = {'max_ambig':0}
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        
    def test_check_ambig_count_detailed(self):
        """ check_ambig_count properly flags for ambiguous base count """
        
        fasta_seq = "ACKCYATASAT"
        final_log_data = {'max_ambig':0}
        detailed_quality_data = {}
        
        failed, final_log_data =\
         check_ambig_count(fasta_seq,
                           final_log_data,
                           detailed_quality_data,
                           max_ambig = 2,
                           enable_all_checks = True)
                           
        expected_failed = True
        expected_final_log_data = {'max_ambig':1}
        expected_detailed_quality_data = {'exceeds_max_ambig':1}
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_check_sliding_qual_window(self):
        """ check_sliding_qual_window truncates low quality window """
        
        fasta_seq = "ACTCGATASAT"
        qual_seq = [10,10,10,10,10,8,8,8,8,8,8]
        final_log_data = {'low_qual_window_found': 1}
        detailed_quality_data = {}
        
        failed, fasta_seq, qual_seq, final_log_data =\
         check_sliding_qual_window(fasta_seq,
                                   qual_seq,
                                   final_log_data,
                                   detailed_quality_data,
                                   qual_score_window = 5,
                                   min_qual_score = 9,
                                   discard_bad_windows = False,
                                   enable_all_checks = False)
        expected_failed = False
        expected_fasta_seq = "ACT"
        expected_qual_seq = [10, 10, 10]
        expected_final_log_data = {'low_qual_window_found': 2}
        expected_detailed_quality_data = {}
        
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_check_sliding_qual_window_detailed(self):
        """ check_sliding_qual_window truncates low quality window """
        
        fasta_seq = "ACTCGATASAT"
        qual_seq = [10,10,10,10,10,8,8,8,8,8,8]
        final_log_data = {'low_qual_window_found': 1,
         'low_qual_window_discarded':0}
        detailed_quality_data = {'discarded_for_low_qual_window':0}
        
        failed, fasta_seq, qual_seq, final_log_data =\
         check_sliding_qual_window(fasta_seq,
                                   qual_seq,
                                   final_log_data,
                                   detailed_quality_data,
                                   qual_score_window = 5,
                                   min_qual_score = 9,
                                   discard_bad_windows = True,
                                   enable_all_checks = True)
        
        expected_failed = True
        expected_fasta_seq = "ACTCGATASAT"
        expected_qual_seq = [10,10,10,10,10,8,8,8,8,8,8]
        expected_final_log_data = {'low_qual_window_discarded': 1,
         'low_qual_window_found': 2}
        expected_detailed_quality_data = {'low_qual_window_found': 1,
         'discarded_for_low_qual_window': 1}
        
        self.assertEqual(fasta_seq, expected_fasta_seq)
        self.assertEqual(qual_seq, expected_qual_seq)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(detailed_quality_data, expected_detailed_quality_data)
        
    def test_check_average_quality(self):
        """ check_average_quality finds average quality score """
        
        qual_seq = [10,10,10,10,10,8,8,8,8,8,8]
        final_log_data = {'min_ave_qual_score': 1}
        detailed_quality_data = {'below_min_ave_qual_score': 1}
        
        failed, final_log_data = check_average_quality(qual_seq,
                          final_log_data,
                          detailed_quality_data,
                          min_qual_score = 9,
                          enable_all_checks = False)
                          
        expected_failed = True
        expected_final_log_data = {'min_ave_qual_score': 2}
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        
    def test_check_average_quality_detailed(self):
        """ check_average_quality finds average quality score """
        
        qual_seq = [10,10,10,10,10,8,8,8,8,8,8]
        final_log_data = {'min_ave_qual_score': 1}
        detailed_quality_data = {}
        
        failed, final_log_data = check_average_quality(qual_seq,
                          final_log_data,
                          detailed_quality_data,
                          min_qual_score = 9,
                          enable_all_checks = True)
                          
        expected_failed = True
        expected_final_log_data = {'min_ave_qual_score': 2}
        expected_detailed_quality_data = {'below_min_ave_qual_score': 1}
        
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(expected_detailed_quality_data, detailed_quality_data)
        
    def test_check_homopolymers(self):
        """ check_homopolymers detects homopolymers in read """
        
        fasta_seq = "ATATTTTCGATTACG"
        final_log_data = {'max_homopolymer': 0}
        detailed_quality_data = {}
        
        failed, final_log_data = check_homopolymers(fasta_seq,
                       final_log_data,
                       detailed_quality_data,
                       max_homopolymer = 3,
                       enable_all_checks = False)
                       
        expected_failed = True
        expected_final_log_data = {'max_homopolymer': 1}
        expected_detailed_quality_data = {}
        
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(expected_detailed_quality_data, detailed_quality_data)
        
    def test_check_homopolymers_detailed(self):
        """ check_homopolymers detects homopolymers in read """
        
        fasta_seq = "ATATTTTCGATTACG"
        final_log_data = {'max_homopolymer': 0}
        detailed_quality_data = {}
        
        failed, final_log_data = check_homopolymers(fasta_seq,
                       final_log_data,
                       detailed_quality_data,
                       max_homopolymer = 3,
                       enable_all_checks = True)
                       
        expected_failed = True
        expected_final_log_data = {'max_homopolymer': 1}
        expected_detailed_quality_data = {'exceeds_max_homopolymer':1}
        
        self.assertEqual(failed, expected_failed)
        self.assertEqual(final_log_data, expected_final_log_data)
        self.assertEqual(expected_detailed_quality_data, detailed_quality_data)
        
    def test_check_map(self):
        """ Should return valid mapping data normally """
        
        header, mapping_data = check_map(self.valid_mapping_data_golay,
              suppress_primer_check = False,
              reverse_primers = 'disable',
              suppress_sampleid_check = False)
              
        expected_header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'ReversePrimer', 'Description']
        expected_mapping_data =\
         [['PC.634', 'AACTCGTCGATG', 'CATGCTTC', 'GACGCGTGAGTAA', 's1_description'],
          ['PC.481', 'agcAGCACTTGT', 'CATGYTGC', 'GTGAGCRACCTGC', 's2_description']]
        
        self.assertEqual(header, expected_header)
        self.assertEqual(mapping_data, expected_mapping_data)
        
    def test_check_map_errors(self):
        """ check_map raises errors only in correct situations """
        
        header_errors_mapping_data =\
         ['#SampleID\tLinkerPrimerSequence\tReversePrimer\tDescription',
         'PC.634\tAACTCGTCGATG\tCATGCTTC\tGACGCGTGAGTAA\ts1_description',
         'PC.633\tagcAGCACTTGT\tCATGYTGC\tGTGAGCRACCTGC\ts2_description']
         
        self.assertRaises(ValueError, check_map, header_errors_mapping_data)
         
        sample_id_errors_mapping_data =\
         ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription',
         'PC.634\tAACTCGTCGATG\tCATGCTTC\tGACGCGTGAGTAA\ts1_description',
         'PC.634\tagcAGCACTTGT\tCATGYTGC\tGTGAGCRACCTGC\ts2_description']
        
        self.assertRaises(ValueError, check_map, sample_id_errors_mapping_data)
        
        # Should not error with a call suppressing the SampleID checks
        header, mapping_data = check_map(sample_id_errors_mapping_data,
              suppress_primer_check = False,
              reverse_primers = 'disable',
              suppress_sampleid_check = True)
              
        primer_errors_mapping_data =\
         ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription',
         'PC.634\tAACTCGTCGATG\tCATXCTTC\tGACGCGTGAGTAA\ts1_description',
         'PC.635\tagcAGCACTTGT\tCATGYTGC\tGTGAGCRACCTGC\ts2_description']
         
        self.assertRaises(ValueError, check_map, primer_errors_mapping_data)
        
        # Should not error with a call suppressing the primer checks
        header, mapping_data = check_map(primer_errors_mapping_data,
              suppress_primer_check = True,
              reverse_primers = 'disable',
              suppress_sampleid_check = False)
        
    def test_get_ids_primers(self):
        """ get_ids_primers generates dictionary of sampleids to primers """
        
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
         'ReversePrimer', 'Description']
        mapping_data = [['PC.634', 'AACTCGTCGATG', 'CATGCTTC',
         'GACGCGTGAGTAA', 's1_description'], ['PC.481', 'agcAGCACTTGT',
         'CATGYTGC', 'GTGAGCRACCTGC', 's2_description']]
        
        ids_primers, ids_rev_primers =\
                    get_ids_primers(header,
                    mapping_data,
                    reverse_primers = "disable",
                    suppress_sampleid_check = False)
                    
        expected_ids_primers = {'PC.634': set(['CATGCTTC']),
         'PC.481': set(['CATGTTGC', 'CATGCTGC'])}
        expected_ids_rev_primers = {}
        
        self.assertEqual(ids_primers, expected_ids_primers)
        self.assertEqual(ids_rev_primers, expected_ids_rev_primers)
        
    def test_get_ids_primers_rev_primers(self):
        """ get_ids_primers also handles reverse primers """
        
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
         'ReversePrimer', 'Description']
        mapping_data = [['PC.634', 'AACTCGTCGATG', 'CATGCTTC',
         'GACGCGTGAGTAA', 's1_description'], ['PC.481', 'agcAGCACTTGT',
         'CATGYTGC', 'GTGAGCRACCTGC', 's2_description']]
        
        ids_primers, ids_rev_primers =\
                    get_ids_primers(header,
                    mapping_data,
                    reverse_primers = "truncate_only",
                    suppress_sampleid_check = False)
                    
        expected_ids_primers = {'PC.634': set(['CATGCTTC']),
         'PC.481': set(['CATGTTGC', 'CATGCTGC'])}
        expected_ids_rev_primers = {'PC.634': set(['TTACTCACGCGTC']),
         'PC.481': set(['GCAGGTTGCTCAC', 'GCAGGTCGCTCAC'])}
        
        self.assertEqual(ids_primers, expected_ids_primers)
        self.assertEqual(ids_rev_primers, expected_ids_rev_primers)
        
    def test_get_ids_primers_suppressed_sampleid_check(self):
        """ get_ids_primers handles suppressed sampleid check """
        
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
         'ReversePrimer', 'Description']
        mapping_data = [['PC.634', 'AACTCGTCGATG', 'CATGCTTC',
         'GACGCGTGAGTAA', 's1_description'], ['PC.481', 'agcAGCACTTGT',
         'CATGYTGC', 'GTGAGCRACCTGC', 's2_description']]
        
        ids_primers, ids_rev_primers =\
                    get_ids_primers(header,
                    mapping_data,
                    reverse_primers = "truncate_only",
                    suppress_sampleid_check = True)
                    
        expected_ids_primers = {'all_primers':
         set(['CATGTTGC', 'CATGCTGC', 'CATGCTTC'])}

        expected_ids_rev_primers =\
         {'rev_primers': set(['GCAGGTTGCTCAC', 'GCAGGTCGCTCAC',
         'TTACTCACGCGTC'])}
        
        self.assertEqual(ids_primers, expected_ids_primers)
        self.assertEqual(ids_rev_primers, expected_ids_rev_primers)
     
    def test_format_detailed_data(self):
        """ Formats detailed log data correctly """
        
        detailed_quality_data = {'seq_order':['seq1', 'seq2'],
         'all_tests':['test1', 'test2'], 'seq1':{'test1':1, 'test2':0},
         'seq2':{'test1':2, 'test2':1}}
        
        formatted_detailed_data =\
         format_detailed_data(detailed_quality_data)
         
        expected_detailed_data = """# Each line contains the fasta label followed by a 1 or 0 for the results of each filter/check
Sequence label	test1	test2
seq1	1	0
seq2	2	1
"""
        self.assertEqual(formatted_detailed_data, expected_detailed_data)
        
        
    def test_write_qual_line(self):
        """ Writes qual line scores correctly """
        
        fake_qual_file = FakeOutFile()
        qual_line = [8,8,8,8,8,8,8,2,2,2,2,2]
        
        write_qual_line(fake_qual_file,
                    label_line = "seq1",
                    qual_seq = qual_line,
                    qual_line_size = 10)
                 
        expected_qual_lines = """>seq1\n8 8 8 8 8 8 8 2 2 2\n2 2\n"""   
        self.assertEqual(fake_qual_file.data, expected_qual_lines)
        
    def test_pair_hmm_align_unaligned_seqs(self):
        """ Returns proper local alignment """
        
        alignment =\
         pair_hmm_align_unaligned_seqs(['AACCGGTT', 'TTTAACCGTTTCGAAC'])
        
        expected_alignment = """>seq_0\nAACCGGTT\n>seq_1\nAACCGTTT\n"""
        
        self.assertEqual(alignment, expected_alignment)
        

    def test_expand_degeneracies(self):
        """generate_possibilities should make possible strings"""
        self.assertEqual(expand_degeneracies(['ACG']), ['ACG'])
        self.assertEqual(expand_degeneracies(['RGY']), 
            ['AGT', 'AGC', 'GGT', 'GGC'])
            
    def test_count_mismatches(self):
        """count_mismatches should count mismatches correctly"""
        self.assertEqual(count_mismatches('GG','GG',0), False)
        self.assertEqual(count_mismatches('GGG','AAA',0), True)
        self.assertEqual(count_mismatches('GGG','AAA',3), False)
        
    def test_count_ambig(self):
        """count_ambig should count ambiguous bases in seq"""
        s = 'ACC'
        s2 = 'RNY'
        s3 = 'NA'
        self.assertEqual(count_ambig(s), 0)
        self.assertEqual(count_ambig(s2), 3)
        self.assertEqual(count_ambig(s3), 1)
        self.assertEqual(count_ambig(''), 0)
        
    def test_make_histograms(self):
        """make_histograms should make correct histograms"""
        raw_lengths = [90, 100, 110, 110, 130, 135]
        post_lengths = [130, 135]
        raw_hist, post_hist, bin_edges = \
            make_histograms(raw_lengths, post_lengths)
        expected_raw_hist = [1, 1, 2, 0, 2, 0]
        expected_post_hist = [0, 0, 0, 0, 2, 0]
        expected_bin_edges = array([ 90, 100, 110, 120, 130, 140, 150])
            
        self.assertEqual(raw_hist, expected_raw_hist)
        self.assertEqual(post_hist, expected_post_hist)
        self.assertEqual(bin_edges, expected_bin_edges)
        
    def test_seq_exceeds_homopolymers(self):
        """seq_exceeds_homopolymers returns True if too many homopolymers"""
        self.assertEqual(seq_exceeds_homopolymers('AAACGA',3), False)
        self.assertEqual(seq_exceeds_homopolymers('AAACGA',2), True)
        self.assertEqual(seq_exceeds_homopolymers('AAACGA',1), True)
        self.assertEqual(seq_exceeds_homopolymers('AAACGATTTT',3), True)
        
    def test_check_window_qual_scores(self):
        """check_window_qual_scores returns False, index if window below qual 
        threshold."""
        scores1 = [8,8,8,8,8,8,8,2,2,2,2,2]
        self.assertEqual(check_window_qual_scores(scores1, 5, 5), (True, 5))
        self.assertEqual(check_window_qual_scores(scores1, 10, 5), (False, 2))
        # windowsize larger than qual score list works
        self.assertEqual(check_window_qual_scores(scores1, 100, 5), (False, 0))
        self.assertEqual(check_window_qual_scores([], 5, 1), True)
        #check each base  in its own window
        self.assertEqual(check_window_qual_scores(scores1, 1, 2), (False, 11))
        self.assertEqual(check_window_qual_scores(scores1, 1, 5), (True, 7))
        
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
        
# Large test data strings at the end for better readability

valid_mapping_data_golay =\
         ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription',
         'PC.634\tAACTCGTCGATG\tCATGCTTC\tGACGCGTGAGTAA\ts1_description', # single mismatch 
         'PC.481\tagcAGCACTTGT\tCATGYTGC\tGTGAGCRACCTGC\ts2_description'] # degenerate char
         
mapping_data_aligned_primer =\
         ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription',
         'PC.634\tAACTCGTCGATG\tCATGCTTC\ts1_description', # single mismatch 
         'PC.481\tagcAGCACTTGT\tTGYTGC\ts2_description'] # primer match 2 bp in
         
no_matching_sampleid_mapping =\
         ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription',
         'PC.XXX\tAACTCGTCGATG\tCATGCTTC\tGACGCGTGAGTAA\ts1_description', # single mismatch 
         'PC.YYY\tagcAGCACTTGT\tCATGYTGC\tGTGAGCRACCTGC\ts2_description'] # degenerate char

valid_fasta_file_no_errors = """>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CATGCTGCCTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG"""

# Will have perfect primer match for third sequence
valid_fasta_file_no_errors2 = """>PC.634_6 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>PC.634_8 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CATGCTTCCTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG"""

# PC.481_4 seq has 4 ambiguous bases added
fasta_file_ambig_bases = """>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CATGCTGCCTCCCGTAGGAGTCNTGGGCYCGTGTCTKCAGTSCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CATGCTGCCTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG"""

fasta_file_diff_lengths = """>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CATGCTGCCTCCCGTAGGAGTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CATGCTGCCTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG"""


valid_qual_file_no_errors = """>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
37 37 37 37 37 37 38 37 33 33 21 21 21 26 33 37 36 36 40 33 24 24 29 33 33 39 39 39 40 39 39 39 40 37 37 37 37 37 37 37 37 37 37 37 32 32 20 20 20 20 20 35 35 37 37 37 37 37 37 37
36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 37 37 37 37 37 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 36 33 28 28 28 28 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 37 28 28 28 37 37 37 37 37 36 36 36 36 36 28 26 26 26 26 28 36 36 36 36 36 36 36 37 38 38 38 38 38 37 37 37 37 37 31 31 31 31 31 31 31
31 31 31 31 31 30 22 22 22 25 25 31 31 31 31 31 31 31 25 25 25 25 25 28
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
37 37 37 37 37 37 40 40 40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36
36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36 36 36 36 36 36 31 31 31 31 31 31 31
>PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
35 35 35 35 35 35 35 35 35 35 23 20 20 31 31 33 33 33 35 23 17 17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35
35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35 35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35
35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33 27 27 25 25 25 27 14 14 14 14 14 25 25 34 34 35 35 35 32 33 33 32 35 35 32 25 25 15 20 20 20 28 35 33 33 33 33 35 35
35 35 35 35 35 35 35 35 35 35 35 35 35 29 24 24 24 29 35 35 35 35 33 33 31 31 34 34 34 34 34 34 31 20 20 20 20 20 31 34 31 31 31 31 32 31 31 33 34 25 25 20 20 18 25 28 28 22 20 22
28 28 28 30 30 29 29 29 30 25 25 25 29 29 26 26 25 23"""

valid_qual_file_no_errors2 = """>PC.634_6 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
37 37 37 37 37 37 38 37 33 33 21 21 21 26 33 37 36 36 40 33 24 24 29 33 33 39 39 39 40 39 39 39 40 37 37 37 37 37 37 37 37 37 37 37 32 32 20 20 20 20 20 35 35 37 37 37 37 37 37 37
36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 37 37 37 37 37 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 36 33 28 28 28 28 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 37 28 28 28 37 37 37 37 37 36 36 36 36 36 28 26 26 26 26 28 36 36 36 36 36 36 36 37 38 38 38 38 38 37 37 37 37 37 31 31 31 31 31 31 31
31 31 31 31 31 30 22 22 22 25 25 31 31 31 31 31 31 31 25 25 25 25 25 28
>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
37 37 37 37 37 37 40 40 40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36
36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36 36 36 36 36 36 31 31 31 31 31 31 31
>PC.634_8 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
35 35 35 35 35 35 35 35 35 35 23 20 20 31 31 33 33 33 35 23 17 17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35
35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35 35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35
35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33 27 27 25 25 25 27 14 14 14 14 14 25 25 34 34 35 35 35 32 33 33 32 35 35 32 25 25 15 20 20 20 28 35 33 33 33 33 35 35
35 35 35 35 35 35 35 35 35 35 35 35 35 29 24 24 24 29 35 35 35 35 33 33 31 31 34 34 34 34 34 34 31 20 20 20 20 20 31 34 31 31 31 31 32 31 31 33 34 25 25 20 20 18 25 28 28 22 20 22
28 28 28 30 30 29 29 29 30 25 25 25 29 29 26 26 25 23"""

# **************** Expected output strings

valid_fasta_demultiplexed = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n']
log_default_demultiplexed = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t1\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/232.00/232.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

valid_fasta_twofiles_demultiplexed = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n', '>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n', '>PC.634_8 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG\n']
log_default_twofile_demultiplexed = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t6\n', 'Seqs written:\t3\n', 'Percent of seqs written\t0.50\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/250.00/238.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t3\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']        
log_min_seq_length_no_seqs_written = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t233\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t0\n', 'Percent of seqs written\t0.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths: NA\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 233:\t1\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

fasta_max_seq_length_changed = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n', '>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n']
log_max_seq_length_changed = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t240\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t6\n', 'Seqs written:\t2\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/232.00/232.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 240:\t1\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t3\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

log_ave_qual_changed = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t35\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t6\n', 'Seqs written:\t2\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/232.00/232.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 35:\t1\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t3\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

fasta_primer_retained = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n']
log_primer_retained = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tTrue\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t1\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t240.00/240.00/240.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

output_fasta_ambig_bases = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCNTGGGCYCGTGTCTKCAGTSCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n']
log_ambig_bases = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t3\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t0\n', 'Percent of seqs written\t0.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths: NA\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 3 ambiguous bases:\t1\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

log_suppress_ambig_check = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t3\n', 'Suppress ambiguous base check:\tTrue\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t1\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/232.00/232.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 3 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

log_changed_homopolymer = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t3\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t0\n', 'Percent of seqs written\t0.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths: NA\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 3:\t1\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

log_suppress_homopolymer = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t3\n', 'Suppress homopolymer check:\tTrue\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t1\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/232.00/232.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 3:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

fasta_output_one_mm_allowed = ['>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA\n', '>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n', '>PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG\n']

log_one_mm_allowed = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t1\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t3\n', 'Percent of seqs written\t1.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/256.00/246.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 1:\t0\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

sliding_window_fasta_out = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCG\n', '>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCG\n']
sliding_window_log = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t50\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t32\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t10\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t6\n', 'Seqs written:\t2\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t162.00/162.00/162.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 50:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 32:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t3\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t1\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t3\n', 'Sequences with detected ambiguous base:\t0\n', '\n']
sliding_window_discard_log = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t50\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t32\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t10\n', 'Discard sequences with low quality windows:\tTrue\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t6\n', 'Seqs written:\t0\n', 'Percent of seqs written\t0.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths: NA\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 50:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 32:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t3\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t3\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t3\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

fasta_out_rev_primer_trunc = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAG\n', '>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAG\n', '>PC.634_8 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG\n']
log_out_rev_primer_trunc = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\ttruncate_only\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t6\n', 'Seqs written:\t3\n', 'Percent of seqs written\t0.50\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t218.00/250.00/228.67\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t3\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t2\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

fasta_out_rev_primer_remove = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAG\n', '>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAG\n']
log_out_rev_primer_remove = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\ttruncate_remove\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t6\n', 'Seqs written:\t2\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t218.00/218.00/218.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t3\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t1\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t2\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

fasta_out_rev_primer_mm = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAG\n', '>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAG\n', '>PC.634_8 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCG\n']
log_out_rev_primer_mm = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\ttruncate_only\n', 'Reverse primer maximum mismatches:\t1\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t6\n', 'Seqs written:\t3\n', 'Percent of seqs written\t0.50\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t218.00/235.00/223.67\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t3\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t3\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

log_out_record_qual_default = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tTrue\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t1\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/232.00/232.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']
qual_out_record_default = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', '40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36 36 36 36 36 36 31 31 31 31 31 31 31\n']

log_out_trunc_ambi = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t10\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tTrue\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tTrue\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t1\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t14.00/14.00/14.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 10:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t1\n', '\n']
qual_out_trunc_rev_primer = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', '40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36\n', '>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', '40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36\n', '>PC.634_8 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', '35 35 23 20 20 31 31 33 33 33 35 23 17 17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26\n', '26 35 35 35 35 35 35 35 33 31 31 31 35 35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33\n', '35 35 35 35 35 35 35 35 33 33 33 35 33 27 27 25 25 25 27 14 14 14 14 14 25 25 34 34 35 35 35 32 33 33 32 35 35 32 25 25 15 20 20 20 28 35 33 33 33 33 35 35 35 35 35 35 35 35 35 35\n', '35 35 35 35 35 29 24 24 24 29 35 35 35 35 33 33 31 31 34 34 34 34 34 34 31 20 20 20 20 20 31 34 31 31 31 31 32 31 31 33 34 25 25 20 20 18 25 28 28 22 20 22 28 28 28 30 30 29 29 29\n', '30 25 25 25 29 29 26 26 25 23\n']
qual_out_trunc_sliding_window = ['>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', '40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '>PC.481_7 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', '40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n', '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37\n']

no_matching_sampleids_fasta = ['>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA\n', '>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n', '>PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CATGCTGCCTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG\n']
no_matching_sampleids_log = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t3\n', 'Percent of seqs written\t1.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t240.00/264.00/254.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t0\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n', 'Fasta labels not matching SampleIDs in mapping file:\n', 'PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n']
sampleids_suppressed_fasta = ['>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA\n', '>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n', '>PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG\n']
sampleids_suppressed_log = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tTrue\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t3\n', 'Percent of seqs written\t1.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/256.00/246.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t0\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

enable_all_tests_main_log = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tTrue\n', 'Record matching quality scores:\tTrue\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t1\n', 'Percent of seqs written\t0.33\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t232.00/232.00/232.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']
enable_all_tests_detailed_log = ['# Each line contains the fasta label followed by a 1 or 0 for the results of each filter/check\n', 'Sequence label\tbelow_min_seq_len\texceeds_max_seq_len\texceeds_max_primer_mismatch\trev_primer_found\texceeds_max_rev_primer_mismatch\tseqs_discarded_for_no_rev_primer\ttoo_short_after_revprimer_truncation\ttruncation_for_ambig_base\ttoo_short_after_ambig_truncation\texceeds_max_ambig\texceeds_max_homopolymer\tbelow_min_ave_qual_score\tlow_qual_window_found\tdiscarded_for_low_qual_window\ttoo_short_after_window_truncation\tseq_id_not_in_mapping\n', 'PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\t0\t0\t1\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\tNA\tNA\tNA\t0\n', 'PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\t0\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\tNA\tNA\tNA\t0\n', 'PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\t0\t0\t1\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\tNA\tNA\tNA\t0\n']

enable_all_tests_main_log_optional_tests = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t35\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t3\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tTrue\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tFalse\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t10\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\ttruncate_only\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tTrue\n', 'Record matching quality scores:\tTrue\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t0\n', 'Percent of seqs written\t0.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths: NA\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 35:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 3 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t2\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t3\n', 'Sequences too short after ambiguous base truncation:\t1\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t1\n', 'Sequences with detected low quality window:\t3\n', 'Sequences with detected ambiguous base:\t1\n', '\n']
enable_all_tests_detailed_log_optional_tests = ['# Each line contains the fasta label followed by a 1 or 0 for the results of each filter/check\n', 'Sequence label\tbelow_min_seq_len\texceeds_max_seq_len\texceeds_max_primer_mismatch\trev_primer_found\texceeds_max_rev_primer_mismatch\tseqs_discarded_for_no_rev_primer\ttoo_short_after_revprimer_truncation\ttruncation_for_ambig_base\ttoo_short_after_ambig_truncation\texceeds_max_ambig\texceeds_max_homopolymer\tbelow_min_ave_qual_score\tlow_qual_window_found\tdiscarded_for_low_qual_window\ttoo_short_after_window_truncation\tseq_id_not_in_mapping\n', 'PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\t0\t0\t1\t0\t1\tNA\t0\t0\tNA\t0\t0\t0\t1\tNA\t1\t0\n', 'PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\t0\t0\t0\t1\t0\tNA\t0\t0\tNA\t0\t0\t0\t1\tNA\t1\t0\n', 'PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\t0\t0\t1\t0\t1\tNA\t0\t1\t1\t0\t0\t0\t1\tNA\t1\t0\n']

expected_fasta_primer_tests_suppressed = ['>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA\n', '>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n', 'CATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG\n', '>PC.634_5 AAAAA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n', 'CATGCTGCCTCCCGTAGGAGTTTGHGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG\n']
expected_log_primer_tests_suppressed = ['\n', 'Parameter Settings\n', '\n', 'Minimum sequence lengths:\t200\n', 'Maximum sequence lengths:\t1000\n', 'Minimum quality score:\t25\n', 'Retain primer sequence:\tFalse\n', 'Maximum ambiguous bases allowed:\t6\n', 'Suppress ambiguous base check:\tFalse\n', 'Truncate at first ambiguous base:\tFalse\n', 'Maximum homopolymer length allowed:\t6\n', 'Suppress homopolymer check:\tFalse\n', 'Suppress primer check:\tTrue\n', 'Maximum allowed primer mismatches:\t0\n', 'Sliding quality window check setting:\t0\n', 'Discard sequences with low quality windows:\tFalse\n', 'Reverse primer settings:\tdisable\n', 'Reverse primer maximum mismatches:\t0\n', 'Suppress SampleID matching:\tFalse\n', 'Enabled detailed logging of all filters:\tFalse\n', 'Record matching quality scores:\tFalse\n', '\n', 'Input sequence count:\t3\n', 'Seqs written:\t3\n', 'Percent of seqs written\t1.00\n', 'Raw min/max/mean sequence lengths:\t240.00/264.00/254.00\n', 'Processed min/max/mean sequence lengths:\t240.00/264.00/254.00\n', '\n', 'Filters that resulted in discarded reads\n', 'Sequences below minimum length 200:\t0\n', 'Sequences above maximum length 1000:\t0\n', 'Sequences below minimum average quality score 25:\t0\n', 'Sequences with homopolymers longer than 6:\t0\n', 'Sequences exceeding 6 ambiguous bases:\t0\n', 'Exceeded maximum primer mismatch of 0:\t0\n', 'Sequences too short after reverse primer truncation:\t0\n', 'Sequences too short after sliding quality window truncation:\t0\n', 'Sequences too short after ambiguous base truncation:\t0\n', 'Sequences discarded for low quality window detection:\t0\n', 'Sequences discarded for failing to detect reverse primer:\t0\n', '\n', 'Other sequence details:\n', '\n', 'Sequences with detected reverse primer:\t0\n', 'Sequences with detected low quality window:\t0\n', 'Sequences with detected ambiguous base:\t0\n', '\n']

if __name__ =='__main__':
    main()