#!/usr/bin/env python
# File created Feb 1 2012

from __future__ import division

__author__ = "William Anton Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Anton Walters"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "William Anton Walters"
__email__ = "william.a.walters@gmail.com"
__status__ = "Development"

from unittest import TestCase, main
from shutil import rmtree
from os.path import exists, join, split

from cogent.util.misc import remove_files

from qiime.util import get_tmp_filename, create_dir
from qiime.validate_demultiplexed_fasta import check_fasta_seqs,\
 get_dup_labels_perc, check_labels_sampleids,\
 run_fasta_checks, validate_fasta, check_fasta_seqs_lens, check_all_ids,\
 check_tree_subset, check_tree_exact_match


class ValidateDemultiplexedFastaTests(TestCase):
    """ """
    
    def setUp(self):
        """ Creates variables and tmp filepaths for use in unit testing """
        
        self.sample_fasta_fp = get_tmp_filename(prefix="sample_fasta_",
         suffix=".fna")
        seq_file = open(self.sample_fasta_fp, 'w')
        seq_file.write(sample_fasta_file)
        seq_file.close()
        
        self.sample_fasta_invalid_fp = get_tmp_filename(prefix="sample_fasta_",
         suffix=".fna")
        seq_file = open(self.sample_fasta_invalid_fp, 'w')
        seq_file.write(sample_fasta_file_invalid)
        seq_file.close()

        self.sample_mapping_fp = get_tmp_filename(prefix="sample_mapping_",
         suffix=".txt")
        map_file = open(self.sample_mapping_fp, "w")
        map_file.write(sample_mapping_file)
        map_file.close()
        
        self.sample_tree_3tips_fp = get_tmp_filename(prefix="sample_tree3tips_",
         suffix=".tre")
        tree_file = open(self.sample_tree_3tips_fp, "w")
        tree_file.write(sample_tree_file_3tips)
        tree_file.close()
        
        self.sample_tree_5tips_fp = get_tmp_filename(prefix="sample_tree3tips_",
         suffix=".tre")
        tree_file = open(self.sample_tree_5tips_fp, "w")
        tree_file.write(sample_tree_file_5tips)
        tree_file.close()
        
        self.sample_mapping_file_errors_fp =\
         get_tmp_filename(prefix="error_mapping_", suffix = ".txt")
        map_file = open(self.sample_mapping_file_errors_fp, "w")
        map_file.write(sample_mapping_file_errors)
        map_file.close() 
         
        self._files_to_remove = [self.sample_fasta_fp,
         self.sample_fasta_invalid_fp, self.sample_mapping_fp,
         self.sample_tree_3tips_fp, self.sample_tree_5tips_fp,
         self.sample_mapping_file_errors_fp]
        
        self.output_dir =\
         get_tmp_filename(prefix = "validate_demultiplexed_fasta_",
         suffix = "/")
        create_dir(self.output_dir)
        
    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)
            
    def test_check_fasta_seqs_lens(self):
        """ Properly returns, sorts sequence lengths """
        
        actual_seq_lens = check_fasta_seqs_lens(self.sample_fasta_fp)
        
        # All seq lengths different
        
        expected_seq_lens = [(1, 35), (1, 32), (1, 27)]
        
        self.assertEqual(actual_seq_lens, expected_seq_lens)
        
        actual_seq_lens = check_fasta_seqs_lens(self.sample_fasta_invalid_fp)
        
        # Should have two seq lens of 35
        
        expected_seq_lens = [(2, 35), (1, 33), (1, 25)]
        
        self.assertEqual(actual_seq_lens, expected_seq_lens)
        
    def test_check_all_ids(self):
        """ Checks that all supplied sampleIDs are in fasta labels """
        
        fasta_labels = ['sample1_1', 'sample1_2', 'sample3_3', 'sample2_4']
        
        sample_ids = ['sample1', 'sample2', 'sample3']
        
        sample_ids_not_found = check_all_ids(fasta_labels, sample_ids)
        
        # should return True as all are found
        
        self.assertEqual(sample_ids_not_found, True)
        
        fasta_labels = ['sample1_1', 'sample1_2', 'sample3_3', 'sample2_4']
        
        sample_ids = ['sample1', 'sample2', 'sample3', 'sampleX']
        
        sample_ids_not_found = check_all_ids(fasta_labels, sample_ids)
        
        # sampleX should not be found
        
        self.assertEqual(sample_ids_not_found, ['sampleX'])
        
        
    def test_check_tree_subset(self):
        """ Checks that fasta labels are a subset of tree tips properly """
        
        fasta_labels = ['seq1_1', 'seq1_2', 'seq2_3', 'seq3_4']
        
        actual_subset_results = check_tree_subset(fasta_labels,
         self.sample_tree_3tips_fp)
         
        # Should find all and give True result
        
        self.assertEqual(actual_subset_results, True)
        
        # Should also get same results with 5 tip tree
        
        fasta_labels = ['seq1_1', 'seq1_2', 'seq2_3', 'seq3_4']
        
        actual_subset_results = check_tree_subset(fasta_labels,
         self.sample_tree_5tips_fp)
         
        # Should find all and give True result
        
        self.assertEqual(actual_subset_results, True)
        
        # Change two of the fasta labels to not match tree tips
        
        fasta_labels = ['seq1_1', 'seqX_2', 'seq2_3', 'seqY_4']
        
        actual_subset_results = check_tree_subset(fasta_labels,
         self.sample_tree_5tips_fp)
         
        # Should find seqX and seqY as not being a subset
        
        self.assertEqual(actual_subset_results, ['seqX', 'seqY'])
        
    def test_check_tree_exact_match(self):
        """ Checks for mismatches between tree tips and fasta labels """
        
        fasta_labels = ['seq1_1', 'seq1_2', 'seq2_3', 'seq3_4']
        
        actual_subset_results = check_tree_exact_match(fasta_labels,
         self.sample_tree_3tips_fp)
         
        # Should find all and give True, True result
        
        self.assertEqual(actual_subset_results, [True, True])
        
        # Should get tips not found in fasta labels with 5 tip tree
        
        fasta_labels = ['seq1_1', 'seq1_2', 'seq2_3', 'seq3_4']
        
        actual_subset_results = check_tree_exact_match(fasta_labels,
         self.sample_tree_5tips_fp)
         
        # Should find all and give True result
        
        self.assertEqual(actual_subset_results, [True, ['seq5', 'seq4']])
        
        # Change two of the fasta labels to not match tree tips
        
        fasta_labels = ['seq1_1', 'seqX_2', 'seq2_3', 'seqY_4']
        
        actual_subset_results = check_tree_exact_match(fasta_labels,
         self.sample_tree_5tips_fp)
         
        # Should find seqX and seqY as not being a subset
        
        self.assertEqual(actual_subset_results, [['seqX', 'seqY'],
         ['seq3', 'seq5', 'seq4']])

    def test_check_fasta_seqs_all_valid(self):
        """ Properly returns invalid chars, barcodes, primers detected """
        
        # Test against all valid data
        
        sample_barcodes = set(['ACCATACC', 'CCAGATTACG'])
        sample_primers = set(['ACATTATTTT', 'TTATTACCGAT'])
        total_seq_count = 3
        
        perc_invalid_chars, perc_barcodes_detected, perc_primers_detected,\
         perc_bcs_seq_start =\
         check_fasta_seqs(self.sample_fasta_fp, sample_barcodes,
         sample_primers, total_seq_count)
         
        expected_perc_invalid_chars = "%1.3f" % 0
        expected_perc_barcodes_detected = "%1.3f" % 0
        expected_perc_primers_detected = "%1.3f" % 0
         
        self.assertEqual(perc_invalid_chars, expected_perc_invalid_chars)
        self.assertEqual(perc_barcodes_detected,
         expected_perc_barcodes_detected)
        self.assertEqual(perc_primers_detected,
         expected_perc_primers_detected)
         
    def test_check_fasta_seqs_with_invalid(self):
        """ Properly returns invalid chars, barcodes, primers detected """
        
        # Test against all data that should give some percent failures
        
        sample_barcodes = set(['ACCATACC', 'AGATTATAT'])
        sample_primers = set(['AGATTTACCA', 'TTATTACCGAT'])
        total_seq_count = 4
        
        perc_invalid_chars, perc_barcodes_detected, perc_primers_detected,\
         perc_bcs_seq_start =\
         check_fasta_seqs(self.sample_fasta_invalid_fp, sample_barcodes,
         sample_primers, total_seq_count)
         
        expected_perc_invalid_chars = "%1.3f" % 0.50
        expected_perc_barcodes_detected = "%1.3f" % 0.25
        expected_perc_primers_detected = "%1.3f" % 0.25
         
        self.assertEqual(perc_invalid_chars, expected_perc_invalid_chars)
        self.assertEqual(perc_barcodes_detected,
         expected_perc_barcodes_detected)
        self.assertEqual(perc_primers_detected,
         expected_perc_primers_detected)
         
    def test_check_fasta_seqs_with_barcodes_at_start(self):
        """ Properly detects barcodes at the start of the sequences """
        
        # Test against all data that should give some percent failures
        
        sample_barcodes = set(['ACAAG', 'AGATTATAT'])
        sample_primers = set(['AGATTTACCA', 'TTATTACCGAT'])
        total_seq_count = 4
        
        perc_invalid_chars, perc_barcodes_detected, perc_primers_detected,\
         perc_bcs_seq_start =\
         check_fasta_seqs(self.sample_fasta_invalid_fp, sample_barcodes,
         sample_primers, total_seq_count)
         
        expected_perc_invalid_chars = "%1.3f" % 0.50
        expected_perc_barcodes_detected = "%1.3f" % 0.50
        expected_perc_primers_detected = "%1.3f" % 0.25
        expected_perc_bcs_seq_start = "%1.3f" % 0.25
         
        self.assertEqual(perc_invalid_chars, expected_perc_invalid_chars)
        self.assertEqual(perc_barcodes_detected,
         expected_perc_barcodes_detected)
        self.assertEqual(perc_primers_detected,
         expected_perc_primers_detected)
        self.assertEqual(perc_bcs_seq_start, expected_perc_bcs_seq_start)
        
        
    def test_get_dup_labels_perc_all_valid(self):
        """ Properly gets percentage of duplicate labels """
        
        # No duplicates
        
        labels = ['seq1', 'seq2', 'seq3', 'seq4']
        
        actual_perc, dups = get_dup_labels_perc(labels)
        
        expected_perc = "%1.3f" % 0.0
        
        self.assertEqual(actual_perc, expected_perc)
        
        expected_dups = []
        
        self.assertEqual(dups, expected_dups)
        
    def test_get_dup_labels_perc_with_dupes(self):
        """ Properly gets percentage of duplicate labels """
        
        # half of the labels are duplicates
        
        labels = ['seq1', 'seq2', 'seq1', 'seq2']
        
        actual_perc, dups = get_dup_labels_perc(labels)
        
        expected_perc = "%1.3f" % 0.5
        
        self.assertEqual(actual_perc, expected_perc)
        
        expected_dups = ['seq1', 'seq2']
        
        self.assertEqual(set(dups), set(expected_dups))
        
    def test_check_labels_sampleids_all_valid(self):
        """  Properly checks for QIIME labels, sampleID matching """
        
        fasta_labels = ['seq1_1', 'seq1_2', 'seq2_3', 'seq1_4', 'seq4_5']
        
        sample_ids = ['seq1', 'seq2', 'seq3', 'seq4']
        
        total_seq_count = 5
        
        perc_not_valid, perc_nosampleid_match =\
         check_labels_sampleids(fasta_labels, sample_ids, total_seq_count)
         
        expected_perc_not_valid = "%1.3f" % 0.0
        expected_perc_nosampleid_match = "%1.3f" % 0.0
        
        self.assertEqual(perc_not_valid, expected_perc_not_valid)
        self.assertEqual(perc_nosampleid_match, expected_perc_nosampleid_match)
        
    def test_check_labels_sampleids_with_problems(self):
        """  Properly checks for QIIME labels, sampleID matching """
        
        fasta_labels = ['seq1_1', 'seq1', 'seq2_3', 'seq1_4', 'seq4_5',
         'seq4_6_1']
        
        sample_ids = ['seq0', 'seq2', 'seq3', 'seq4']
        
        total_seq_count = 6
        
        perc_not_valid, perc_nosampleid_match =\
         check_labels_sampleids(fasta_labels, sample_ids, total_seq_count)
         
        expected_perc_not_valid = "%1.3f" % 0.333
        expected_perc_nosampleid_match = "%1.3f" % 0.667
        
        self.assertEqual(perc_not_valid, expected_perc_not_valid)
        self.assertEqual(perc_nosampleid_match, expected_perc_nosampleid_match)
        
    def test_run_fasta_checks(self):
        """ Properly returns dictionary of percentage for each check """
        
        actual_fasta_report = run_fasta_checks(self.sample_fasta_fp,
         self.sample_mapping_fp)
         
        # All values should be zero
        expected_fasta_report = {'invalid_labels': '0.000',
         'tree_subset': False, 
         'all_ids_found': False, 
         'same_seq_lens': False, 
         'barcodes_detected': '0.000', 
         'duplicate_labels': '0.000',
         'duplicate_ids': [],
         'invalid_seq_chars': '0.000', 
         'nosample_ids_map': '0.000', 
         'linkerprimers_detected': '0.000', 
         'tree_exact_match': False,
         'barcodes_at_start': '0.000'}
         
        self.assertEqual(actual_fasta_report, expected_fasta_report)
        
    def test_run_fasta_checks_with_invalid_data(self):
        """ Properly returns dictionary of percentage for each check """
        
        actual_fasta_report = run_fasta_checks(self.sample_fasta_invalid_fp,
         self.sample_mapping_fp)
         
        # All values should be nonzero
        expected_fasta_report = {'invalid_labels': '0.500',
         'tree_subset': False,
         'all_ids_found': False,
         'same_seq_lens': False, 
         'barcodes_detected': '0.250', 
         'duplicate_labels': '0.250',
         'duplicate_ids': ['seq1'], 
         'invalid_seq_chars': '0.500', 
         'nosample_ids_map': '0.750', 
         'linkerprimers_detected': '0.250', 
         'tree_exact_match': False,
         'barcodes_at_start': '0.000'}
         
         
        self.assertEqual(actual_fasta_report, expected_fasta_report)
        
    def test_validate_fasta_suppress_primers_barcodes(self):
        """ Overall module test with primer/barcode check suppressed """
       
        # Should raise errors when both primer/barcode check not suppressed
        self.assertRaises(ValueError, validate_fasta, self.sample_fasta_fp,
         self.sample_mapping_file_errors_fp, self.output_dir)
         
        self.assertRaises(ValueError, validate_fasta, self.sample_fasta_fp,
         self.sample_mapping_file_errors_fp, self.output_dir,
         suppress_primer_checks=True)
        
        self.assertRaises(ValueError, validate_fasta, self.sample_fasta_fp,
         self.sample_mapping_file_errors_fp, self.output_dir,
         suppress_barcode_checks=True)
        
        # No errors when both suppressed
        validate_fasta(self.sample_fasta_fp, self.sample_mapping_file_errors_fp,
         self.output_dir, suppress_primer_checks = True,
         suppress_barcode_checks = True)
         
        expected_log_fp = join(self.output_dir,
         split(self.sample_fasta_fp)[1] + "_report.log")
         
        log_f = open(expected_log_fp, "U")
        actual_log_lines = [line.strip() for line in log_f][1:]
        
        expected_log_lines = """Percent duplicate labels: 0.000
Percent QIIME-incompatible fasta labels: 0.000
Percent of labels that fail to map to SampleIDs: 0.000
Percent of sequences with invalid characters: 0.000
Percent of sequences with barcodes detected: 0.000
Percent of sequences with barcodes detected at the beginning of the sequence: 0.000
Percent of sequences with primers detected: 0.000""".split('\n')

        self.assertEqual(actual_log_lines, expected_log_lines)
        
    def test_validate_fasta(self):
        """ Overall module runs properly """
        
                   
        validate_fasta(self.sample_fasta_fp, self.sample_mapping_fp,
         self.output_dir)
         
        expected_log_fp = join(self.output_dir,
         split(self.sample_fasta_fp)[1] + "_report.log")
         
        log_f = open(expected_log_fp, "U")
        actual_log_lines = [line.strip() for line in log_f][1:]
        
        expected_log_lines = """Percent duplicate labels: 0.000
Percent QIIME-incompatible fasta labels: 0.000
Percent of labels that fail to map to SampleIDs: 0.000
Percent of sequences with invalid characters: 0.000
Percent of sequences with barcodes detected: 0.000
Percent of sequences with barcodes detected at the beginning of the sequence: 0.000
Percent of sequences with primers detected: 0.000""".split('\n')

        self.assertEqual(actual_log_lines, expected_log_lines)
        
        # Check with all optional values included
        
        validate_fasta(self.sample_fasta_fp, self.sample_mapping_fp,
         self.output_dir, tree_fp=self.sample_tree_5tips_fp, tree_subset=True,
         tree_exact_match=True, same_seq_lens=True, all_ids_found=True)
         
        expected_log_fp = join(self.output_dir,
         split(self.sample_fasta_fp)[1] + "_report.log")
         
        log_f = open(expected_log_fp, "U")
        actual_log_lines = [line.strip() for line in log_f][1:]
        
        expected_log_lines = """Percent duplicate labels: 0.000
Percent QIIME-incompatible fasta labels: 0.000
Percent of labels that fail to map to SampleIDs: 0.000
Percent of sequences with invalid characters: 0.000
Percent of sequences with barcodes detected: 0.000
Percent of sequences with barcodes detected at the beginning of the sequence: 0.000
Percent of sequences with primers detected: 0.000
Sequence lengths report
Counts of sequences, followed by their sequence lengths:
1\t35
1\t32
1\t27
Sample ID in fasta sequences report
The following SampleIDs were not found:
seq2
Fasta label subset in tree tips report
All fasta labels were a subset of tree tips.
Fasta label/tree tip exact match report
All fasta labels found in tree tips.
The following tips were not in fasta labels:
seq2
seq5
seq4""".split('\n')

        self.assertEqual(actual_log_lines, expected_log_lines)
        
    def test_validate_fasta_with_invalid(self):
        """ Overall module runs properly """
                   
        validate_fasta(self.sample_fasta_invalid_fp, self.sample_mapping_fp,
         self.output_dir)
         
        expected_log_fp = join(self.output_dir,
         split(self.sample_fasta_invalid_fp)[1] + "_report.log")
         
        log_f = open(expected_log_fp, "U")
        actual_log_lines = [line.strip() for line in log_f][1:]
        
        expected_log_lines = """Percent duplicate labels: 0.250
Percent QIIME-incompatible fasta labels: 0.500
Percent of labels that fail to map to SampleIDs: 0.750
Percent of sequences with invalid characters: 0.500
Percent of sequences with barcodes detected: 0.250
Percent of sequences with barcodes detected at the beginning of the sequence: 0.000
Percent of sequences with primers detected: 0.250
Duplicate labels found:
seq1""".split('\n')

        self.assertEqual(actual_log_lines, expected_log_lines)

# Should all be valid
sample_fasta_file = """>seq1_1
ACAAGCAGAGATAATTACCAGAGTTTT
>seq1_2
ACCAGGAGATTTAACCGGAGAGGAGAGAGTTT
>seq3_3
ACAGAGTTTATATAGGAGAGATTTAACCAGNAACC"""

# some invalid data
sample_fasta_file_invalid =""">seq1
ACAAGCAGAGATTTACCAGAGTTTT
>seq1
ACCAGGAGATTTAACCGGAGAGGAGAGAGTTTX
>seq3_3
ACAGAGTTTATATAGGAGAGATTTAACCAGNAACC
>seq4_5
ACCAGATTATATTACCAGATTACCCAGAG--ATTA"""


# sample mapping file
sample_mapping_file = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription
seq1\tAAAAAAAA\tTTTTTTTT\tseq1_desc
seq2\tTTTTTTTT\tAGATTTACCA\tseq2_desc
seq3\tAGATTATAT\tTTTTTTYT\tseq3_desc"""

sample_mapping_file_errors = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription
seq1\t\t\tseq1_desc
seq2\t\tAGATTTACCA\tseq2_desc
seq3\tAGATTATAT\tTTTTTTYT\tseq3_desc"""

# sample tree file
sample_tree_file_5tips = """(seq1:0.18882,seq2:0.14833,(seq3:0.11088,(seq4:0.1312,seq5:0.11842)0.890:0.01354)0.752:0.0132);"""

sample_tree_file_3tips = """(seq1:0.18882,seq2:0.14833,(seq3:0.11088)0.752:0.0132);"""

if __name__ == "__main__":
    main()
