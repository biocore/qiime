#!/usr/bin/env python

"""
provides unit tests for the usearch.py module
"""

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["William Walters", 
               "Jose Carlos Clemente Litran",
               "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"

from os.path import isfile, basename, join, exists
from shutil import rmtree
from glob import glob

from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationError, get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import create_dir, get_random_directory_name, remove_files
from qiime.pycogent_backports.usearch import (Usearch,
 clusters_from_blast_uc_file, usearch_fasta_sort_from_filepath,
 usearch_dereplicate_exact_subseqs, usearch_dereplicate_exact_seqs,
 usearch_sort_by_abundance,
 usearch_cluster_error_correction, usearch_chimera_filter_de_novo,
 usearch_chimera_filter_ref_based, usearch_cluster_seqs,
 enumerate_otus, assign_reads_to_otus, usearch_qf, concatenate_fastas,
 get_retained_chimeras, assign_dna_reads_to_protein_database,
 assign_dna_reads_to_dna_database, usearch61_ref_cluster,
 usearch61_denovo_cluster, sort_by_abundance_usearch61,
 sort_by_length_usearch61, usearch61_cluster_ref,
 usearch61_fast_cluster, usearch61_smallmem_cluster,
 parse_dereplicated_uc, parse_usearch61_clusters,
 merge_clusters_dereplicated_seqs, merge_failures_dereplicated_seqs,
 parse_usearch61_failures, usearch61_chimera_check_denovo,
 usearch61_chimera_check_ref)


class Usearch61Tests(TestCase):
    """ Tests for usearch 6.1 functionality """
    
    def setUp(self):
        # create the temporary input files
        
        self.output_dir = '/tmp/'
        
        self.dna_seqs_1 = dna_seqs_1
        self.usearch_ref_seqs1 = usearch_ref_seqs1
        self.dna_seqs_1_subset = dna_seqs_1_subset
        self.dna_seqs_with_dups = dna_seqs_with_dups2
        self.usearch61_dereplicated_uc_lines = usearch61_dereplicated_uc_lines
        self.usearch61_clustered_uc_lines = usearch61_clustered_uc_lines
        self.usearch61_clustered_uc_lines_ref =\
         usearch61_clustered_uc_lines_ref
        self.usearch61_clustered_ref_lines = usearch61_clustered_ref_lines
        self.de_novo_chimera_seqs = de_novo_chimera_seqs
        self.expected_usearch61_denovo_uchime_file =\
         expected_usearch61_denovo_uchime_file
        self.reference_seqs_fp = reference_seqs_fp
        self.expected_usearch61_ref_uchime_file =\
         expected_usearch61_ref_uchime_file 
        
        self.tmp_dna_seqs_1 = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_dna_seqs_1,'w')
        seq_file.write(self.dna_seqs_1)
        seq_file.close()
        
        self.tmp_usearch_ref_seqs1 = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_usearch_ref_seqs1,'w')
        seq_file.write(self.usearch_ref_seqs1)
        seq_file.close()
        
        self.tmp_dna_seqs_1_subset = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_dna_seqs_1_subset,'w')
        seq_file.write(self.dna_seqs_1_subset)
        seq_file.close()
        
        self.tmp_dna_seqs_with_dups = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_dna_seqs_with_dups, "w")
        seq_file.write(self.dna_seqs_with_dups)
        seq_file.close()
        
        self.tmp_de_novo_chimera_seqs = get_tmp_filename(\
         prefix='Usearch61denovoChimera_',\
         suffix='.fasta')
        seq_file = open(self.tmp_de_novo_chimera_seqs, 'w')
        seq_file.write(self.de_novo_chimera_seqs)
        seq_file.close()
        
        self.tmp_ref_chimera_seqs = get_tmp_filename(\
         prefix="Usearch61refChimera_",\
         suffix='.fasta')
        seq_file = open(self.tmp_ref_chimera_seqs, "w")
        seq_file.write(self.reference_seqs_fp)
        seq_file.close()
                 
        self._files_to_remove =\
         [self.tmp_dna_seqs_1, self.tmp_usearch_ref_seqs1,
          self.tmp_dna_seqs_1_subset, self.tmp_dna_seqs_with_dups,
          self.tmp_de_novo_chimera_seqs, self.tmp_ref_chimera_seqs]
          
        self._dirs_to_remove = []
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        if self._dirs_to_remove:
            for curr_dir in self._dirs_to_remove:
                rmtree(curr_dir)
                
    def test_usearch61_ref_default_params(self):
        """ usearch61 reference OTU picking works with default settings """
        
        clusters, failures = usearch61_ref_cluster(self.tmp_dna_seqs_1,
         self.tmp_usearch_ref_seqs1, output_dir = self.output_dir,
         save_intermediate_files = False, remove_usearch_logs = True)
         
        # Should all fall into single, de novo clusters
         
        expected_failures = []
        
        self.assertEqual(failures, expected_failures)
        
        expected_clusters = [['uclust_test_seqs_9'], ['uclust_test_seqs_8'],
         ['uclust_test_seqs_3'], ['uclust_test_seqs_5'], ['uclust_test_seqs_4'],
         ['uclust_test_seqs_1'], ['uclust_test_seqs_0'], ['uclust_test_seqs_2'],
         ['uclust_test_seqs_7'], ['uclust_test_seqs_6']]
         
        self.assertEqual(len(expected_clusters), 10)
        
        for curr_cluster in clusters.values():
            self.assertTrue(curr_cluster in expected_clusters)

        
    def test_usearch61_ref_default_params_suppressed_clusters(self):
        """ usearch61 reference OTU picking, suppressed clusters """
        
        clusters, failures = usearch61_ref_cluster(self.tmp_dna_seqs_1,
         self.tmp_usearch_ref_seqs1, suppress_new_clusters=True,
         output_dir = self.output_dir,
         save_intermediate_files = False, remove_usearch_logs = True)
         
        # Should all fail as the reference database does not match.
        
        expected_clusters = {}
         
        expected_failures = ['uclust_test_seqs_0', 'uclust_test_seqs_9',
         'uclust_test_seqs_4', 'uclust_test_seqs_7', 'uclust_test_seqs_2',
         'uclust_test_seqs_1', 'uclust_test_seqs_3', 'uclust_test_seqs_8',
         'uclust_test_seqs_6', 'uclust_test_seqs_5']
        
        self.assertEqual(clusters, expected_clusters)
        
        for curr_failure in failures:
            self.assertTrue(curr_failure in expected_failures)
            
    def test_usearch61_ref_default_params_matches_ref(self):
        """ usearch61 reference OTU picking, matches ref OTU IDs """
        
        clusters, failures = usearch61_ref_cluster(self.tmp_dna_seqs_1,
         self.tmp_dna_seqs_1, suppress_new_clusters=True,
         output_dir = self.output_dir,
         save_intermediate_files = False, remove_usearch_logs = True)
         
        # Should all fall into single, ref-based clusters
        
        expected_clusters = {'uclust_test_seqs_5': ['uclust_test_seqs_5'],
         'uclust_test_seqs_4': ['uclust_test_seqs_4'],
         'uclust_test_seqs_7': ['uclust_test_seqs_7'],
         'uclust_test_seqs_6': ['uclust_test_seqs_6'],
         'uclust_test_seqs_1': ['uclust_test_seqs_1'],
         'uclust_test_seqs_0': ['uclust_test_seqs_0'],
         'uclust_test_seqs_3': ['uclust_test_seqs_3'],
         'uclust_test_seqs_2': ['uclust_test_seqs_2'],
         'uclust_test_seqs_9': ['uclust_test_seqs_9'],
         'uclust_test_seqs_8': ['uclust_test_seqs_8']}
         
        expected_failures = []
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
        
    def test_usearch61_ref_open_ref(self):
        """ usearch61 does open reference OTU picking """
        
        clusters, failures = usearch61_ref_cluster(self.tmp_dna_seqs_1,
         self.tmp_dna_seqs_1_subset, percent_id = 0.98, rev = True,
         save_intermediate_files = False, minlen = 44,
         output_dir = self.output_dir, remove_usearch_logs = True,
         verbose = False, wordlength = 12, usearch_fast_cluster = False,
         usearch61_sort_method = 'abundance', otu_prefix = "denovo",
         usearch61_maxrejects = 100, usearch61_maxaccepts = 4,
         sizeorder=True)
         
        # Should all fall into single, ref-based & denovo clusters
        
        expected_ref_results = {'uclust_test_seqs_1': ['uclust_test_seqs_1'],
         'uclust_test_seqs_0': ['uclust_test_seqs_0'],
         'uclust_test_seqs_3': ['uclust_test_seqs_3'],
         'uclust_test_seqs_2': ['uclust_test_seqs_2']}
         
        expected_denovo_results = [['uclust_test_seqs_5'], 
         ['uclust_test_seqs_7'], ['uclust_test_seqs_8'], ['uclust_test_seqs_4'],
         ['uclust_test_seqs_6'], ['uclust_test_seqs_9']]
         
        self.assertEqual(len(clusters), 10)
        
        for curr_ref_result in expected_ref_results:
            self.assertEqual(clusters[curr_ref_result],
             expected_ref_results[curr_ref_result])
        for curr_denovo_result in expected_denovo_results:
            self.assertTrue(curr_denovo_result in clusters.values())
         
        expected_failures = []

        self.assertEqual(failures, expected_failures)
        
    def test_usearch61_denovo_default_params(self):
        """ usearch61 denovo OTU picking works with default settings """
        
        clusters = usearch61_denovo_cluster(self.tmp_dna_seqs_1,
         output_dir = self.output_dir, save_intermediate_files = False,
         remove_usearch_logs = True)
         
        # Should all fall into single, de novo clusters
        
        expected_clusters = [['uclust_test_seqs_9'], ['uclust_test_seqs_8'],
         ['uclust_test_seqs_3'], ['uclust_test_seqs_5'], ['uclust_test_seqs_4'],
         ['uclust_test_seqs_1'], ['uclust_test_seqs_0'], ['uclust_test_seqs_2'],
         ['uclust_test_seqs_7'], ['uclust_test_seqs_6']]
         
        self.assertEqual(len(expected_clusters), 10)
        
        for curr_cluster in clusters.values():
            self.assertTrue(curr_cluster in expected_clusters)
        
    def test_usearch61_denovo_length_sorting(self):
        """ usearch61 denovo OTU picking works with length sorting """
        
        clusters = usearch61_denovo_cluster(self.tmp_dna_seqs_1,
         output_dir = self.output_dir, save_intermediate_files = False,
         remove_usearch_logs = True, usearch61_sort_method = 'length')
         
        # Should all fall into single, de novo clusters
        
        expected_clusters = [['uclust_test_seqs_9'], ['uclust_test_seqs_8'],
         ['uclust_test_seqs_3'], ['uclust_test_seqs_5'], ['uclust_test_seqs_4'],
         ['uclust_test_seqs_1'], ['uclust_test_seqs_0'], ['uclust_test_seqs_2'],
         ['uclust_test_seqs_7'], ['uclust_test_seqs_6']]
         
        self.assertEqual(len(expected_clusters), 10)
        
        for curr_cluster in clusters.values():
            self.assertTrue(curr_cluster in expected_clusters)
        
    def test_usearch61_denovo_no_sorting(self):
        """ usearch61 denovo OTU picking works with no sorting """
        
        clusters = usearch61_denovo_cluster(self.tmp_dna_seqs_1,
         output_dir = self.output_dir, save_intermediate_files = False,
         remove_usearch_logs = True, usearch61_sort_method = 'None')
         
        # Should all fall into single, de novo clusters
        
        expected_clusters = [['uclust_test_seqs_9'], ['uclust_test_seqs_8'],
         ['uclust_test_seqs_3'], ['uclust_test_seqs_5'], ['uclust_test_seqs_4'],
         ['uclust_test_seqs_1'], ['uclust_test_seqs_0'], ['uclust_test_seqs_2'],
         ['uclust_test_seqs_7'], ['uclust_test_seqs_6']]
         
        self.assertEqual(len(expected_clusters), 10)
        
        for curr_cluster in clusters.values():
            self.assertTrue(curr_cluster in expected_clusters)
        
    def test_usearch61_denovo_fast_cluster(self):
        """ usearch61 denovo OTU picking works with fast_cluster sorting """
        
        clusters = usearch61_denovo_cluster(self.tmp_dna_seqs_1,
         output_dir = self.output_dir, save_intermediate_files = False,
         remove_usearch_logs = True, usearch61_sort_method = 'length',
         usearch_fast_cluster = True)
         
        # Should all fall into single, de novo clusters
        
        expected_clusters = [['uclust_test_seqs_9'], ['uclust_test_seqs_8'],
         ['uclust_test_seqs_3'], ['uclust_test_seqs_5'], ['uclust_test_seqs_4'],
         ['uclust_test_seqs_1'], ['uclust_test_seqs_0'], ['uclust_test_seqs_2'],
         ['uclust_test_seqs_7'], ['uclust_test_seqs_6']]
         
        self.assertEqual(len(expected_clusters), 10)
        
        for curr_cluster in clusters.values():
            self.assertTrue(curr_cluster in expected_clusters)
            
    def test_sort_by_abundance_usearch61(self):
        """ usearch61 sorts by abundance successfully """
        
        sorted_fna_fp = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        sorted_uc_fp = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.uc')
        
        output_fna_filepath, output_uc_filepath, app_result =\
         sort_by_abundance_usearch61(self.tmp_dna_seqs_with_dups,
         self.output_dir, remove_usearch_logs = True,
         output_fna_filepath = sorted_fna_fp,
         output_uc_filepath = sorted_uc_fp, log_name = "abundance_sorted.log")
         
        output_fna = [line for line in MinimalFastaParser(open(output_fna_filepath, "U"))]
        
        expected_fna = [('seq2;size=3;', 
         'TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC'),
         ('seq1;size=1;',
         'GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAA')]        
        
        self._files_to_remove.append(sorted_fna_fp)
        self._files_to_remove.append(sorted_uc_fp)
        
        self.assertEqual(output_fna, expected_fna)
        
    def test_sort_by_length_usearch61(self):
        """ usearch61 sorts by length successfully """
        
        sorted_fna_fp = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        
        output_fna_filepath, app_result =\
         sort_by_length_usearch61(self.tmp_usearch_ref_seqs1,
         self.output_dir, remove_usearch_logs = True,
         output_fna_filepath = sorted_fna_fp)
         
        output_fna = [line for line in MinimalFastaParser(open(output_fna_filepath, "U"))]
        
        expected_fna = [('ref1',
         'CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCA'),
         ('L07864',
         'GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTAATGCATGGGAATCTGCCATATAGTGGGGGACAACTGGGGAAACCCAGGCTAATACCGCATAATCTCTACGGAGGAAAGGCTTC'),
         ('EU199232',
         'TACGCGCGGAAATCGAGCGAGATTGGGAACGCAAGTTCCTGAGTATTGCGGCGAACGGGTGAGTAAGACGTGGGTGATCTACCCCTAGGGTGGGAATAACCCGGGGAAACCCGGGCTAATACCGAATAAGACCACAGGAGGCGACTCCAGAGGGTCAAAGGGAGCCTTGGCCTCCCCC')]        
        self._files_to_remove.append(sorted_fna_fp)
        
        self.assertEqual(output_fna, expected_fna)
        
    def test_usearch61_cluster_ref(self):
        """ usearch61 reference OTU picking application call successful """
        
        output_uc_fp = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.uc')
        
        uc_fp, failures = usearch61_cluster_ref(self.tmp_dna_seqs_1,
         self.tmp_dna_seqs_1, output_dir = self.output_dir,
         remove_usearch_logs = True, output_uc_filepath=output_uc_fp)
         
        self._files_to_remove.append(uc_fp)
         
        actual_uc_lines = [line.strip() for line in open(uc_fp, "U")]
        
        # Difficult to test output, as numbers change between runs, for now
        # just testing length of output, and order of lines changes as well.
        
        self.assertEqual(len(actual_uc_lines), 10)
        
    def test_usearch61_fast_cluster(self):
        """ usearch61 fast cluster OTU picking application call successful """
        
        output_uc_fp = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.uc')
        
        uc_fp, failures = usearch61_fast_cluster(self.tmp_dna_seqs_1,
         output_dir = self.output_dir,
         remove_usearch_logs = True, output_uc_filepath=output_uc_fp)
         
        self._files_to_remove.append(uc_fp)
         
        actual_uc_lines = [line.strip() for line in open(uc_fp, "U")]
        
        # Difficult to test output, as numbers change between runs, for now
        # just testing length of output, and order of lines changes as well.
        
        self.assertEqual(len(actual_uc_lines), 20)
        
    def test_usearch61_cluster_smallmem(self):
        """ usearch61 smallmem OTU picking application call successful """
        
        output_uc_fp = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.uc')
        
        uc_fp, failures = usearch61_smallmem_cluster(self.tmp_dna_seqs_1,
         output_dir = self.output_dir,
         remove_usearch_logs = True, output_uc_filepath=output_uc_fp)
         
        self._files_to_remove.append(uc_fp)
         
        actual_uc_lines = [line.strip() for line in open(uc_fp, "U")]
        
        # Difficult to test output, as numbers change between runs, for now
        # just testing length of output, and order of lines changes as well.
        
        self.assertEqual(len(actual_uc_lines), 20)
        
    def test_parse_dereplicated_uc(self):
        """ Parses dereplicated usearch61 uc file successfully """
        
        actual_derep_ids =\
         parse_dereplicated_uc(self.usearch61_dereplicated_uc_lines)

        expected_derep_ids = {'seq2': ['seq3', 'seq4'], 'seq1': []}
         
        self.assertEqual(actual_derep_ids, expected_derep_ids)
        
    def test_parse_usearch61_clusters_denovo(self):
        """ Parses usearch61 de novo clusters uc file correctly """
        
        actual_parsed_clusters, failures =\
         parse_usearch61_clusters(self.usearch61_clustered_uc_lines,
         ref_clustered=False)
        
         
        expected_parsed_clusters =\
         ({'denovo0': ['seq2'], 'denovo1': ['seq1']})
         
        self.assertEqual(actual_parsed_clusters, expected_parsed_clusters)
        
    def test_parse_usearch61_clusters_ref(self):
        """ Parses usearch61 ref clusters uc file correctly """
        
        actual_parsed_clusters, failures =\
         parse_usearch61_clusters(self.usearch61_clustered_uc_lines_ref,
         otu_prefix = '', ref_clustered=True)
         
        expected_parsed_clusters =\
         ({'seq4': ['seq2'], 'seq1': ['seq1']})
         
        self.assertEqual(actual_parsed_clusters, expected_parsed_clusters)
        
    def test_merge_clusters_dereplicated_seqs(self):
        """ Properly merges dereplicated and clustered sequences """
        
        derep_ids = {'seq2': ['seq3', 'seq4'], 'seq1': []}
        
        clustered_ids = ({'seq4': ['seq2'], 'seq1': ['seq1']})
        
        merged_ids = merge_clusters_dereplicated_seqs(clustered_ids,
         derep_ids)
         
        expected_ids = {'seq1': ['seq1'], 'seq4': ['seq2', 'seq3', 'seq4']}
        
        self.assertEqual(merged_ids, expected_ids)
        
    def test_merge_failures_dereplicated_seqs(self):
        """ Usearch61 properly merges dereplicated seqs, ref based failures """
        
        failures = ['seq2']
        derep_ids = {'seq2': ['seq3', 'seq4'], 'seq1': []}
        
        merged_failures = merge_failures_dereplicated_seqs(failures,
         derep_ids)
         
        expected_failures = ['seq2', 'seq3', 'seq4']
        
        self.assertEqual(merged_failures, expected_failures)
        
    def test_parse_usearch61_failures(self):
        """ Writes failures out to fasta file """
        
        failures = ['seq2', 'seq3', 'seq4']
        filtered_fna_fp = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        output_fp = parse_usearch61_failures(self.tmp_dna_seqs_with_dups,
         failures, filtered_fna_fp)
         
        self._files_to_remove.append(output_fp)
         
        output_fna = [line for line in MinimalFastaParser(open(output_fp, "U"))]

        expected_fna = [('seq2', 'TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC'), ('seq3', 'TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC'), ('seq4', 'TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC')]
        self.assertEqual(output_fna, expected_fna)
        
    # Chimera tests
    
    def test_usearch61_denovo_chimera_detection(self):
        """ usearch61 denovo chimera detection correctly flags chimeras """
        
        uchime_fp = join(self.output_dir, "uchime_denovo.uchime")
        
        
        uchime_fp, app_result =\
         usearch61_chimera_check_denovo(self.tmp_de_novo_chimera_seqs,
                                        uchime_denovo_fp = uchime_fp,
                                        output_dir = self.output_dir,
                                        remove_usearch_logs = True)
                                        
        uchime_f = open(uchime_fp, "U")
        
        actual_lines = [line.strip() for line in uchime_f]
        
        # There is some system dependent stochastic effect on calculations
        # for chimeras, need to pull out only the flags Y or N for chimeras
          
        expected_chimera_ixs = [11, 16]
        
        for line in range(len(actual_lines)):
            curr_chimera_flag = actual_lines[line].split('\t')[-1]
            if line in expected_chimera_ixs:
                self.assertEqual(curr_chimera_flag, "Y")
            else:
                self.assertEqual(curr_chimera_flag, "N")
                        

        self._files_to_remove.append(uchime_fp)
        

        
    def test_usearch61_ref_chimera_detection(self):
        """ usearch61 ref chimera detection correctly flags chimeras """
        
        uchime_fp = join(self.output_dir, "uchime_ref.uchime")
        
        
        uchime_fp, app_result =\
         usearch61_chimera_check_ref(self.tmp_de_novo_chimera_seqs,
                                     uchime_ref_fp = uchime_fp,
                                     reference_seqs_fp =\
                                      self.tmp_ref_chimera_seqs,
                                     output_dir = self.output_dir,
                                     remove_usearch_logs = True)
                                        
        uchime_f = open(uchime_fp, "U")
        
        actual_lines = [line.strip() for line in uchime_f]
        
        self.assertEqual(actual_lines,
         self.expected_usearch61_ref_uchime_file)

        self._files_to_remove.append(uchime_fp)
        
        
        
        
class UsearchTests(TestCase):

    def setUp(self):
        # create the temporary input files
        self.dna_seqs_1 = dna_seqs_1
        self.dna_seqs_2 = dna_seqs_usearch
        self.dna_seqs_3 = dna_seqs_3
        self.dna_seqs_4 = dna_seqs_4
        self.protein_ref_seqs1 = protein_ref_seqs1
        self.ref_database = usearch_ref_seqs1
        self.dna_seqs_with_abundance = dna_seqs_with_abundance
        self.de_novo_chimera_seqs = de_novo_chimera_seqs
        self.dna_seqs_with_dups = dna_seqs_with_dups
        self.dna_seqs_reference_otu_picking = dna_seqs_reference_otu_picking
        
        
        # Expected output files
        self.uc_lines1 = uc_lines1
        self.expected_otu_assignments = expected_otu_assignments
        self.expected_enumerated_fasta = expected_enumerated_fasta
        self.expected_enumerated_fasta_added_options =\
         expected_enumerated_fasta_added_options
        self.expected_clusters_w_abundance_default_settings =\
         expected_clusters_w_abundance_default_settings
        self.expected_clusters_w_abundance_low_setting =\
         expected_clusters_w_abundance_low_setting
        self.expected_reference_filtered_seqs =\
         expected_reference_filtered_seqs
        self.expected_de_novo_chimeras_default =\
         expected_de_novo_chimeras_default
        self.expected_de_novo_chimera_filtered_skew11 =\
         expected_de_novo_chimera_filtered_skew11
        self.expected_cluster_err_seqs =\
         expected_cluster_err_seqs
        self.expected_sorted_by_abundance_no_filter =\
         expected_sorted_by_abundance_no_filter
        self.expected_derep_seqs = expected_derep_seqs
        self.expected_abundance_sort_filtered = expected_abundance_sort_filtered
        self.expected_len_sorted_seqs = expected_len_sorted_seqs
        self.expected_combined_dna_seqs_1_seqs_usearch =\
         expected_combined_dna_seqs_1_seqs_usearch
        self.retained_chimeras_seqs1 = retained_chimeras_seqs1
        self.retained_chimeras_seqs2 = retained_chimeras_seqs2
        self.expected_retained_chimeras_union =\
         expected_retained_chimeras_union
        self.expected_retained_chimeras_intersection =\
         expected_retained_chimeras_intersection
        self.expected_derep_seqs_full_len =\
         expected_derep_seqs_full_len
        
        # Create temporary files for use with unit tests
        
        self.tmp_dir = '/tmp/'
        
        self.tmp_seq_filepath1 = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1,'w')
        seq_file.write(self.dna_seqs_1)
        seq_file.close()        
        
        self.tmp_seq_filepath2 = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath2,'w')
        seq_file.write(self.dna_seqs_2)
        seq_file.close()

        self.dna_seqs3_filepath = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.dna_seqs3_filepath,'w')
        seq_file.write(self.dna_seqs_3)
        seq_file.close()

        self.dna_seqs4_filepath = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.dna_seqs4_filepath,'w')
        seq_file.write(self.dna_seqs_4)
        seq_file.close()

        self.protein_ref_seqs1_filepath = get_tmp_filename(\
         prefix='UsearchOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.protein_ref_seqs1_filepath,'w')
        seq_file.write(self.protein_ref_seqs1)
        seq_file.close()
        
        self.tmp_ref_database = get_tmp_filename(\
         prefix='UsearchRefDatabase_',\
         suffix='.fasta')
        seq_file = open(self.tmp_ref_database, 'w')
        seq_file.write(self.ref_database)
        seq_file.close()
        
        self.tmp_seqs_w_abundance = get_tmp_filename(\
         prefix='UsearchSeqsAbundance_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seqs_w_abundance, 'w')
        seq_file.write(self.dna_seqs_with_abundance)
        seq_file.close()
        
        self.tmp_de_novo_chimera_seqs = get_tmp_filename(\
         prefix='UsearchdenovoChimera_',\
         suffix='.fasta')
        seq_file = open(self.tmp_de_novo_chimera_seqs, 'w')
        seq_file.write(self.de_novo_chimera_seqs)
        seq_file.close()
        
        self.tmp_dna_seqs_with_dups = get_tmp_filename(\
         prefix='UsearchDupDNASeqs_',\
         suffix='.fasta')
        seq_file = open(self.tmp_dna_seqs_with_dups, 'w')
        seq_file.write(self.dna_seqs_with_dups)
        seq_file.close()
        
        self.tmp_retained_chimeras_seqs1 = get_tmp_filename(\
         prefix="UsearchRetainedChimeras1_",\
         suffix=".fasta")
        seq_file = open(self.tmp_retained_chimeras_seqs1, 'w')
        seq_file.write(self.retained_chimeras_seqs1)
        seq_file.close()
        
        self.tmp_retained_chimeras_seqs2 = get_tmp_filename(\
         prefix="UsearchRetainedChimeras1_",\
         suffix=".fasta")
        seq_file = open(self.tmp_retained_chimeras_seqs2, 'w')
        seq_file.write(self.retained_chimeras_seqs2)
        seq_file.close()
        
        self.tmp_dna_seqs_ref_otu_picking = get_tmp_filename(\
         prefix="UsearchRefOtuPicking_",\
         suffix=".fasta")
        seq_file = open(self.tmp_dna_seqs_ref_otu_picking, "w")
        seq_file.write(self.dna_seqs_reference_otu_picking)
        seq_file.close()
        
       
        self._files_to_remove =\
         [self.tmp_seq_filepath1, self.tmp_seq_filepath2,
          self.tmp_ref_database, self.tmp_seqs_w_abundance,
          self.tmp_de_novo_chimera_seqs, self.tmp_dna_seqs_with_dups,
          self.tmp_retained_chimeras_seqs1, self.tmp_retained_chimeras_seqs2,
          self.tmp_dna_seqs_ref_otu_picking, self.dna_seqs3_filepath,
          self.protein_ref_seqs1_filepath, self.dna_seqs4_filepath]
          
        self._dirs_to_remove = []
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        if self._dirs_to_remove:
            for curr_dir in self._dirs_to_remove:
                rmtree(curr_dir)
            
    def test_usearch_qf(self):
        """ Main program loop test, with default parameters """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = usearch_qf(self.tmp_seq_filepath2,
                                      output_dir = self.tmp_dir,
                                      db_filepath = self.tmp_ref_database,
                                      minsize = 1,
                                      remove_usearch_logs=True,
                                      chimeras_retention = 'intersection')
                                      
        expected_clusters = {'1': ['Solemya', 'Solemya_seq2'],
         '0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2']}
        expected_failures = ['chimera']
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
    
    def test_usearch_qf_minlen(self):
        """ Main program loop test, with longer minlen """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = usearch_qf(self.tmp_seq_filepath2,
                                      output_dir = self.tmp_dir,
                                      db_filepath = self.tmp_ref_database,
                                      minsize = 1,
                                      remove_usearch_logs=True,
                                      chimeras_retention = 'intersection',
                                      minlen=110)
                                      
        expected_clusters = {'0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2']}
        expected_failures = ['Solemya', 'Solemya_seq2', 'chimera']
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
        
    def test_usearch_qf_reference_otu_picking(self):
        """ Main program loop test, with reference + new clusters """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = usearch_qf(self.tmp_dna_seqs_ref_otu_picking,
                                      output_dir = self.tmp_dir,
                                      refseqs_fp = self.tmp_ref_database,
                                      reference_chimera_detection=False,
                                      minsize = 1,
                                      remove_usearch_logs=True,
                                      suppress_new_clusters=False)
                                      
        # Will cluster everything including RandomCrap, as new clusters allowed.
        expected_clusters =  {'1': ['Solemya', 'Solemya_seq2'], 
        '0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2'],
        '2': ['RandomCrap']}
        expected_failures = []
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
    
    def test_usearch_qf_reference_otu_picking_no_new_clusters(self):
        """ Main program loop test, with reference and no new clusters """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = usearch_qf(self.tmp_dna_seqs_ref_otu_picking,
                                      output_dir = self.tmp_dir,
                                      refseqs_fp = self.tmp_ref_database,
                                      reference_chimera_detection=False,
                                      minsize = 1,
                                      remove_usearch_logs=True,
                                      suppress_new_clusters=True)
                                      
        # Will cluster everything but RandomCrap, as no new clusters allowed.
        expected_clusters =  {'L07864': ['Solemya', 'Solemya_seq2'], 
        'ref1': ['usearch_ecoli_seq', 'usearch_ecoli_seq2']}
        expected_failures = ['RandomCrap']
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
        
        
    def test_usearch_qf_no_ref_database(self):
        """ Main program loop with no reference chimera testing """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = usearch_qf(self.tmp_seq_filepath2,
                                      output_dir = self.tmp_dir,
                                      reference_chimera_detection=False,
                                      minsize = 1,
                                      remove_usearch_logs=True)
                                      
        # Chimera sequence should not be detected without reference test.
        expected_clusters = {'1': ['Solemya', 'Solemya_seq2'],
         '0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2'],
         '2': ['chimera']}
         
        expected_failures = []
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
        
    def test_usearch_qf_union(self):
        """ Main program loop with union nonchimera retention """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = usearch_qf(self.tmp_seq_filepath2,
                                      output_dir = self.tmp_dir,
                                      reference_chimera_detection=False,
                                      minsize = 1,
                                      remove_usearch_logs=True,
                                      chimeras_retention = 'union')
                                      
        # Chimera sequence retained as passes de novo test
        expected_clusters = {'1': ['Solemya', 'Solemya_seq2'],
         '0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2'],
         '2': ['chimera']}
         
        expected_failures = []
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)

        
    def test_usearch_qf_disabled_filters(self):
        """ Returns expected clustering with no filtering """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = usearch_qf(self.tmp_seq_filepath2,
                                      output_dir = self.tmp_dir,
                                      de_novo_chimera_detection=False,
                                      reference_chimera_detection=False,
                                      cluster_size_filtering=False,
                                      remove_usearch_logs=True)
                                      
        # Chimera sequence should not be detected without reference test.
        expected_clusters = {'1': ['Solemya', 'Solemya_seq2'],
         '0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2'],
         '2': ['chimera']}
         
        expected_failures = []
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
        
    def test_usearch_qf_generates_logs(self):
        """ Generates expected log files """
        
        curr_output_dir = get_tmp_filename(prefix='/UsearchLogTest_',suffix='/')
         
        create_dir(curr_output_dir)
        
         
        self._dirs_to_remove.append(curr_output_dir)
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = usearch_qf(self.tmp_seq_filepath2,
                                      output_dir = curr_output_dir,
                                      db_filepath = self.tmp_ref_database,
                                      minsize = 1,
                                      remove_usearch_logs=False,
                                      chimeras_retention = 'intersection')
                                      
        expected_clusters = {'1': ['Solemya', 'Solemya_seq2'],
         '0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2']}
        expected_failures = ['chimera']
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
        
        # Only checking for creation of files, as file contents contain 
        # tmp file names.
        expected_log_names = ['assign_reads_to_otus.log',
                              'uchime_de_novo_chimera_filtering.log',
                              'derep.log',
                              'uchime_reference_chimera_filtering.log',
                              'minsize_0_abundance_sort.log',
                              'usearch_cluster_err_corrected.log',
                              'minsize_1_abundance_sort.log',
                              'usearch_cluster_seqs.log',
                              'sortlen.log']
                              
        actual_logs =\
         [basename(curr_file) for curr_file in glob(curr_output_dir + "*.*")]
        
        for log in expected_log_names:
            self.assertContains(actual_logs, log)
            
    def test_concatenate_fastas(self):
        """ Properly concatenates two fasta files """
        
        out_f =\
         get_tmp_filename(prefix='UsearchConcatFileTest_',suffix='.fasta')
        
        actual_concatenated_seqs = concatenate_fastas(self.tmp_seq_filepath1,
         self.tmp_seq_filepath2, out_f)
        
        self._files_to_remove.append(out_f)
        
        actual_lines =\
         [line.strip() for line in open(actual_concatenated_seqs, "U")]
         
        self.assertEqual(actual_lines,
         expected_combined_dna_seqs_1_seqs_usearch)
         
        
    def test_assign_reads_to_otus(self):
        """ Properly assigns reads back to original ID """
        
        
        app_result, output_filepath =\
         assign_reads_to_otus(original_fasta = self.tmp_ref_database,
                              filtered_fasta = self.tmp_seq_filepath2,
                              remove_usearch_logs = True,
                              working_dir=self.tmp_dir)
                              
        self._files_to_remove.append(output_filepath)
        
        # Stripping off first line, which refers to the command using tmp
        # file names, retaining other actual results.
        actual_assignments =\
         [line.strip() for line in open(output_filepath, "U")][2:]
        
        self.assertEqual(actual_assignments, self.expected_otu_assignments)
                              
        
    def test_enumerate_otus(self):
        """ Enumerates OTUs properly """
        
        output_filepath = enumerate_otus(self.tmp_seq_filepath1)
        
        self._files_to_remove.append(output_filepath)

        actual_fasta = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_fasta, self.expected_enumerated_fasta)
        
    def test_enumerate_otus_added_options(self):
        """ Enumerates with all options properly """
        
        output_filepath = enumerate_otus(self.tmp_seq_filepath1,
                                         label_prefix = "Big",
                                         label_suffix = "Ern",
                                         retain_label_as_comment = True,
                                         count_start = 255)
        
        self._files_to_remove.append(output_filepath)

        actual_fasta = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_fasta,
         self.expected_enumerated_fasta_added_options)
        
        
    def test_usearch_cluster_seqs(self):
        """ Clusters sequences correctly """
        
        # clusters all seqs with default 97% identity
        app_result, output_filepath =\
         usearch_cluster_seqs(self.tmp_seqs_w_abundance,
                              save_intermediate_files=False,
                              remove_usearch_logs=True,
                              percent_id = 0.97,
                              working_dir=self.tmp_dir)
        
        self._files_to_remove.append(output_filepath)
                              
        actual_clusters = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_clusters,
         self.expected_clusters_w_abundance_default_settings)
        
                            
    def test_usearch_cluster_seqs_high_identity(self):
        """ Clusters sequences correctly """
        
        # Should get two clusters with 99.9% identity
        app_result, output_filepath =\
         usearch_cluster_seqs(self.tmp_seqs_w_abundance,
                              save_intermediate_files=False,
                              remove_usearch_logs=True,
                              percent_id = 0.999,
                              working_dir=self.tmp_dir)
        
        self._files_to_remove.append(output_filepath)
                              
        actual_clusters = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_clusters,
         self.expected_clusters_w_abundance_low_setting)
        
    def test_usearch_chimera_filter_ref_based(self):
        """ Properly detects chimeras against reference database """
        
        app_result, output_filepath =\
         usearch_chimera_filter_ref_based(self.tmp_seq_filepath2,
                                          self.tmp_ref_database,
                                          remove_usearch_logs=True,
                                          working_dir=self.tmp_dir)
                                          
        self._files_to_remove.append(output_filepath)
        
        actual_filtered_chimeras =\
         [line.strip() for line in open(output_filepath, "U")]
         
        self.assertEqual(actual_filtered_chimeras,
         self.expected_reference_filtered_seqs)
                                          
        
    def test_usearch_chimera_filter_de_novo(self):
        """ Properly detects de novo chimeras """
        
        app_result, output_filepath =\
         usearch_chimera_filter_de_novo(self.tmp_de_novo_chimera_seqs,
                                        remove_usearch_logs=True,
                                        abundance_skew = 2,
                                        working_dir=self.tmp_dir)
                                        
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = \
         [line.strip() for line in open(output_filepath, "U")]
         
        self.assertEqual(actual_seqs, self.expected_de_novo_chimeras_default)
        
    def test_usearch_chimera_filter_de_novo_abundance_skew(self):
        """ Properly detects de novo chimeras with skew changes """
        
        app_result, output_filepath =\
         usearch_chimera_filter_de_novo(self.tmp_de_novo_chimera_seqs,
                                        remove_usearch_logs=True,
                                        abundance_skew = 11,
                                        working_dir=self.tmp_dir)
                                        
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = \
         [line.strip() for line in open(output_filepath, "U")]
         
        self.assertEqual(actual_seqs,
         self.expected_de_novo_chimera_filtered_skew11)
        
    def test_usearch_cluster_error_correction(self):
        """ Properly clusters seqs for chimera testing/filtering """
        
        # clusters all seqs with default 97% identity
        app_result, output_filepath =\
         usearch_cluster_error_correction(self.tmp_seqs_w_abundance,
                                          save_intermediate_files=False,
                                          remove_usearch_logs=True,
                                          percent_id_err = 0.97,
                                          working_dir=self.tmp_dir)
        
        self._files_to_remove.append(output_filepath)
                              
        actual_clusters = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_clusters,
         self.expected_cluster_err_seqs)
        
    def test_usearch_sort_by_abundance(self):
        """ Properly sorts fasta by abundance """
        
        app_result, output_filepath =\
         usearch_sort_by_abundance(self.tmp_de_novo_chimera_seqs,
                                   remove_usearch_logs=True,
                                   working_dir=self.tmp_dir)
                                   
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_seqs,
         self.expected_sorted_by_abundance_no_filter)
         
    def test_usearch_sort_by_abundance_filter(self):
        """ Properly sorts fasta by abundance, filters low count otus """
        
        app_result, output_filepath =\
         usearch_sort_by_abundance(self.tmp_de_novo_chimera_seqs,
                                   remove_usearch_logs=True,
                                   minsize = 40,
                                   working_dir=self.tmp_dir)
                                   
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_seqs,
         self.expected_abundance_sort_filtered)

    def test_usearch_dereplicate_exact_subseqs(self):
        """ Properly dereplicates fasta file """
        
        app_result, output_filepath =\
         usearch_dereplicate_exact_subseqs(self.tmp_dna_seqs_with_dups,
                                           remove_usearch_logs=True,
                                           working_dir=self.tmp_dir)
                                           
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = [line.strip() for line in open(output_filepath, "U")]

        
        self.assertEqual(actual_seqs, self.expected_derep_seqs)

    def test_usearch_dereplicate_exact_seqs(self):
        """ Properly dereplicates fasta file """

        app_result, output_filepath =\
                    usearch_dereplicate_exact_seqs(self.tmp_dna_seqs_with_dups,
                                                      remove_usearch_logs=True,
                                                      working_dir=self.tmp_dir)

        self._files_to_remove.append(output_filepath)

        actual_seqs = [line.strip() for line in open(output_filepath, "U")]

        self.assertEqual(actual_seqs, self.expected_derep_seqs_full_len)
        
    def test_usearch_fasta_sort_from_filepath(self):
        """ Properly sorts fasta according to seq length """
        
        app_result, output_filepath =\
         usearch_fasta_sort_from_filepath(self.tmp_seq_filepath2,
                                          remove_usearch_logs=True,
                                          working_dir=self.tmp_dir)
                                          
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_seqs, self.expected_len_sorted_seqs)
        

    def test_clusters_from_blast_uc_file(self):
        """ clusters_from_uc_file functions as expected """

        expected_clusters = {'19': ['PC.634_4'], '42': ['PC.test2_1',
         'PC.test1_2', 'PC.634_3'], '6': ['PC.269_5']}
        expected_failures = ['PC.481_6']

        self.assertEqual(clusters_from_blast_uc_file(self.uc_lines1),
         (expected_clusters,expected_failures))
         
    def test_get_retained_chimeras_union(self):
        """ Properly returns union of two fastas """
        
        out_f =\
         get_tmp_filename(prefix='UsearchUnionTest_',suffix='.fasta')
         
        actual_out_fp = get_retained_chimeras(self.tmp_retained_chimeras_seqs1,
         self.tmp_retained_chimeras_seqs2, out_f, chimeras_retention = 'union')
         
        self._files_to_remove.append(out_f)
         
        actual_out_f = [line.strip() for line in open(actual_out_fp, "U")]
        
        self.assertEqual(actual_out_f, self.expected_retained_chimeras_union)
        
    def test_get_retained_chimeras_intersection(self):
        """ Properly returns intersection of two fastas """
        
        out_f =\
         get_tmp_filename(prefix='UsearchIntersectionTest_',suffix='.fasta')
         
        actual_out_fp = get_retained_chimeras(self.tmp_retained_chimeras_seqs1,
         self.tmp_retained_chimeras_seqs2, out_f,
         chimeras_retention = 'intersection')
         
        self._files_to_remove.append(out_f)
         
        actual_out_f = [line.strip() for line in open(actual_out_fp, "U")]
        
        self.assertEqual(actual_out_f,
         self.expected_retained_chimeras_intersection)

    def test_assign_dna_reads_to_protein_database(self):
        """assign_dna_reads_to_protein_database wrapper functions as expected 
        """
        output_dir = get_random_directory_name(output_dir=self.tmp_dir)
        self._dirs_to_remove.append(output_dir)
        output_fp = join(output_dir,'out.uc')
        assign_dna_reads_to_protein_database(self.dna_seqs3_filepath, 
                                         self.protein_ref_seqs1_filepath, 
                                         output_fp, 
                                         temp_dir = self.tmp_dir)
        
        self.assertTrue(exists(output_fp))
        self.assertTrue(exists(output_fp.replace('.uc','.bl6')))
        
        # confirm that the clusters look like what we expect
        expected_clusters = [['eco:b0015'],['eco:b0122','eco:b0122-like']]
        expected_clusters.sort()
        actual_clusters = clusters_from_blast_uc_file(open(output_fp))[0].values()
        actual_clusters.sort()
        self.assertEqual(actual_clusters,expected_clusters)


    def test_assign_dna_reads_to_protein_database_alt_params(self):
        """assign_dna_reads_to_protein_database wrapper functions with alt params 
        """
        output_dir = get_random_directory_name(output_dir=self.tmp_dir)
        self._dirs_to_remove.append(output_dir)
        output_fp = join(output_dir,'out.uc')
        assign_dna_reads_to_protein_database(self.dna_seqs3_filepath, 
                                         self.protein_ref_seqs1_filepath, 
                                         output_fp, 
                                         temp_dir = self.tmp_dir,
                                         params={'--id':1.0})
        
        self.assertTrue(exists(output_fp))
        self.assertTrue(exists(output_fp.replace('.uc','.bl6')))
        
        # confirm that the clusters look like what we expect
        expected_clusters = [['eco:b0015'],['eco:b0122']]
        expected_clusters.sort()
        actual_clusters = clusters_from_blast_uc_file(open(output_fp))[0].values()
        actual_clusters.sort()
        self.assertEqual(actual_clusters,expected_clusters)

    def test_assign_dna_reads_to_dna_database(self):
        """assign_dna_reads_to_protein_database wrapper functions as expected 
        """
        output_dir = get_random_directory_name(output_dir=self.tmp_dir)
        self._dirs_to_remove.append(output_dir)
        output_fp = join(output_dir,'out.uc')
        assign_dna_reads_to_protein_database(self.dna_seqs3_filepath, 
                                         self.dna_seqs4_filepath, 
                                         output_fp, 
                                         temp_dir = self.tmp_dir)
        
        self.assertTrue(exists(output_fp))
        self.assertTrue(exists(output_fp.replace('.uc','.bl6')))
        
        # confirm that the clusters look like what we expect
        expected_clusters = [['eco:b0015'],['eco:b0122','eco:b0122-like']]
        expected_clusters.sort()
        actual_clusters = clusters_from_blast_uc_file(open(output_fp))[0].values()
        actual_clusters.sort()
        self.assertEqual(actual_clusters,expected_clusters)

# Long strings for test files, output, etc.
# *************************************************

retained_chimeras_seqs1 = """>seq1
ACAGGCC
>seq2
ACAGGCCCCC
>seq3
TTATCCATT"""

retained_chimeras_seqs2 = """>seq3
ACAGGCC
>seq4
ACAGGCCCCC
>seq5
TTATCCATT"""

dna_seqs_1 = """>uclust_test_seqs_0 some comment0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_1 some comment1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_2 some comment2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3 some comment3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>uclust_test_seqs_4 some comment4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>uclust_test_seqs_5 some comment4_again
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>uclust_test_seqs_6 some comment6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>uclust_test_seqs_7 some comment7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_8 some comment8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>uclust_test_seqs_9 some comment9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA"""

dna_seqs_1_subset = """>uclust_test_seqs_0 some comment0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_1 some comment1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_2 some comment2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3 some comment3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT"""

dna_seqs_3 = """>eco:b0001 thrL; thr operon leader peptide; K08278 thr operon leader peptide (N)
atgaaacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcg
ggctga
>eco:b0015 dnaJ; chaperone Hsp40, co-chaperone with DnaK; K03686 molecular chaperone DnaJ (N)
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
>eco:b0122 yacC; conserved protein, PulS_OutS family (N)
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa
>eco:b0122-like
atgaagaaaattttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaatcc"""

dna_seqs_4 = """>eco:b0015 dnaJ; chaperone Hsp40, co-chaperone with DnaK; K03686 molecular chaperone DnaJ (N)
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
>eco:b0122 yacC; conserved protein, PulS_OutS family (N)
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa
>eco:b0122-like
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaatcc"""

protein_ref_seqs1 = """>eco:b0001 thrL; thr operon leader peptide; K08278 thr operon leader peptide (A)
MKRISTTITTTITITTGNGAG
>eco:b0015 dnaJ; chaperone Hsp40, co-chaperone with DnaK; K03686 molecular chaperone DnaJ (A)
MAKQDYYEILGVSKTAEEREIRKAYKRLAMKYHPDRNQGDKEAEAKFKEIKEAYEVLTDS
QKRAAYDQYGHAAFEQGGMGGGGFGGGADFSDIFGDVFGDIFGGGRGRQRAARGADLRYN
MELTLEEAVRGVTKEIRIPTLEECDVCHGSGAKPGTQPQTCPTCHGSGQVQMRQGFFAVQ
QTCPHCQGRGTLIKDPCNKCHGHGRVERSKTLSVKIPAGVDTGDRIRLAGEGEAGEHGAP
AGDLYVQVQVKQHPIFEREGNNLYCEVPINFAMAALGGEIEVPTLDGRVKLKVPGETQTG
KLFRMRGKGVKSVRGGAQGDLLCRVVVETPVGLNERQKQLLQELQESFGGPTGEHNSPRS
KSFFDGVKKFFDDLTR
>eco:b0015:rep
MAKQDYYEILGVSKTAEEREIRKAYKRLAMKYHPDRNQGDKEAEAKFKEIKEAYEVLTDS
QKRAAYDQYGHAAFEQGGMGGGGFGGGADFSDIFGDVFGDIFGGGRGRQRAARGADLRYN
MELTLEEAVRGVTKEIRIPTLEECDVCHGSGAKPGTQPQTCPTCHGSGQVQMRQGFFAVQ
QTCPHCQGRGTLIKDPCNKCHGHGRVERSKTLSVKIPAGVDTGDRIRLAGEGEAGEHGAP
AGDLYVQVQVKQHPIFEREGNNLYCEVPINFAMAALGGEIEVPTLDGRVKLKVPGETQTG
KLFRMRGKGVKSVRGGAQGDLLCRVVVETPVGLNERQKQLLQELQESFGGPTGEHNSPRS
KSFFDGVKKFFDDLTR
>eco:b0122 yacC; conserved protein, PulS_OutS family (A)
MKTFFRTVLFGSLMAVCANSYALSESEAEDMADLTAVFVFLKNDCGYQNLPNGQIRRALV
FFAQQNQWDLSNYDTFDMKALGEDSYRDLSGIGIPVAKKCKALARDSLSLLAYVK"""

usearch_ref_seqs1 = """>ref1 ecoli sequence
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCA
>EU199232 1 1236 Bacteria/Deltaproteobacteria/Desulfurella - Hippea/uncultured
TACGCGCGGAAATCGAGCGAGATTGGGAACGCAAGTTCCTGAGTATTGCGGCGAACGGGTGAGTAAGACGTGGGTGATCTACCCCTAGGGTGGGAATAACCCGGGGAAACCCGGGCTAATACCGAATAAGACCACAGGAGGCGACTCCAGAGGGTCAAAGGGAGCCTTGGCCTCCCCC
>L07864 1 1200 Bacteria/Beta Gammaproteobacteria/Solemya symbiont
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTAATGCATGGGAATCTGCCATATAGTGGGGGACAACTGGGGAAACCCAGGCTAATACCGCATAATCTCTACGGAGGAAAGGCTTC
"""

dna_seqs_usearch = """>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGT
>Solemya seq
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTA
>usearch_ecoli_seq2
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTCCAT
>Solemya_seq2
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTATCAAG
>chimera
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACCCCTAGGGTGGGAATAACCCGGGGAAACCCGGGCTAATACCGAATAAGACCACAGGAGGCGACTCCAGAGGGTCAAAGGGAGCCTTGGCCTCCCCC
"""

dna_seqs_reference_otu_picking = """>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGT
>Solemya seq
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTA
>usearch_ecoli_seq2
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTCCAT
>Solemya_seq2
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTATCAAG
>RandomCrap
ACACAAACAGTATATTATATCCCCAGACAGGGACCGAGATTTACCACACCCAAAAAAAAAAAAAACACACCCCCCCCCCCCCCACACACACACTTATTTT
"""

dna_seqs_with_abundance = """>Cluster1;size=114
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Cluster2;size=45
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCC
>Cluster0;size=37
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGAACATCGCGTCAGGTTTGTGTCAGGCCT
>Cluster7;size=33
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAAGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Cluster6;size=32
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Cluster5;size=25
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Cluster11;size=22
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Cluster12;size=15
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Cluster13;size=2
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTGTGTCAGGCCT
>Cluster14;size=1
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTG"""

de_novo_chimera_seqs = """>Cluster1;size=52
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACAC
CGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGT
TATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster0;size=50
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCC
CGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCA
TCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGG
>Cluster2;size=45
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCC
CGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>Cluster10;size=43
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCC
CGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTC
CTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>Cluster4;size=40
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>Cluster6;size=40
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTA
CGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGG
>Cluster3;size=30
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCT
TACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTG
TTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>Cluster12;size=19
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTGGTCTTGGTGAGCCGTTACCC
CACCAACTAACTAATACAACGCATGCCCATCCATTACCACCGGAGTTTTCAACCCAAGAAGATGCCTCCCTGGATGTTAT
GGGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTAATGGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTC
GCCGACAAT
>Cluster30;size=18
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCC
CGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAAGAGATGCCTCTCGGAGTTAT
TACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTACGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster29;size=18
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGACGCCGCGTC
ACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
TCGCCGG
>Cluster16;size=16
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCC
CTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTA
TGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACT
>Cluster222;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAT
GCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACT
>Cluster221;size=1
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTGCCC
CGCCAACTACCTAATCGGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGATTATCTCACCATGCGGCAAAATAATGT
CATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGCGTT
>Cluster218;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGTCACCTTGGTAGGCCGTTACCC
CACCAACTAGCTAATCAGACGCAAGCCCATCTATCAGCGGATTGCTCCTTTTCTAGCTATATCATGCGATACTACTAGCT
TATGCGGTATTAGCAATGATTTCTCACTGTTATTCCCCTCTGATAGGCAGG
>Cluster217;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGATATCCAGAGCCATGCGACCCAGATATATT
ATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCATCCAAGGCAGGTT
>Cluster216;size=1
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGTCAGCTACTGATCGTCGCCTTGGTAGGCCATTACCC
TACCAACTAGCTAATCAGACGCGAGGCCATCTCTCAGCGATAAATCTTTGATATATCTGCCATGCGACAAACATATATTA
TGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTT
>Cluster522;size=10
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAT
GCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACT"""

reference_seqs_fp = """>seq1
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACAC
CGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGT
TATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>seq2
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCC
CGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCA
TCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGG
>seq3
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCC
CGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>mixed_seq
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCAACATATTTCGGGACAGATTAACACACAAAGGATTTACACAAAAT
ACATTAGACCAAACCCCAAGATTTAGACAGGATTACAGGATTTACAGATTTTTACCAACATTAGACAGGGG"""

dna_seqs_with_dups=""">seq1
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAA
>seq2
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
>seq3
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
>seq4
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTT"""

dna_seqs_with_dups2=""">seq1
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAA
>seq2
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
>seq3
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
>seq4
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC"""


# Expected output file data
uc_lines1 = """# usearch --id 0.97 --uc usearch_picked_otus/assign_reads_to_otus.uc --query seqs.fna --global --db usearch_picked_otus/enumerated_otus.fasta
# version=4.2.66
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
H\t42\t217\t99.1\t+\t0\t0\t217MI\tPC.test2_1 FLP3FBN02xELBSXx orig_bc=ACAGAGTCGGCG new_bc=ACAGAGTCGGCG,FLP3FBN02x bc_diffs=0\t42
H\t42\t217\t99.1\t+\t0\t0\t217MI\tPC.test1_2 FLP3FBN03ELBSXx orig_bc=ACAGAGTCGGCG new_bc=ACAGAGTCGGCG,FLP3FBN03 bc_diffs=0\t42
H\t42\t217\t99.1\t+\t0\t0\t217MI\tPC.634_3 FLP3FBN01ELBSX orig_bc=TCAGAGTCGGCT new_bc=ACAGAGTCGGCT,FLP3FBN01 bc_diffs=1\t42
H\t19\t243\t100.0\t+\t0\t0\t25MI218M\tPC.634_4 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT,FLP3FBN01 bc_diffs=0\t19
N\t*\t219\t*\t*\t*\t*\t*\tPC.481_6\tFLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG,FLP3FBN01 bc_diffs=0\t*
H\t6\t211\t99.5\t+\t0\t0\t211M\tPC.269_5 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA,FLP3FBN01 bc_diffs=0\t6
""".split('\n')

expected_otu_assignments = """# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
H\t2\t199\t97.5\t.\t0\t0\t119M80D\tref1 ecoli sequence\tusearch_ecoli_seq2
N\t*\t178\t*\t*\t*\t*\t*\tEU199232 1 1236 Bacteria/Deltaproteobacteria/Desulfurella - Hippea/uncultured\t*
H\t1\t180\t100.0\t.\t0\t0\t97M83D\tL07864 1 1200 Bacteria/Beta Gammaproteobacteria/Solemya symbiont\tSolemya seq""".split('\n')

expected_enumerated_fasta = """>0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>5
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA""".split('\n')

expected_enumerated_fasta_added_options = """>Big255Ern\tuclust_test_seqs_0 some comment0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Big256Ern\tuclust_test_seqs_1 some comment1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>Big257Ern\tuclust_test_seqs_2 some comment2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>Big258Ern\tuclust_test_seqs_3 some comment3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>Big259Ern\tuclust_test_seqs_4 some comment4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>Big260Ern\tuclust_test_seqs_5 some comment4_again
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>Big261Ern\tuclust_test_seqs_6 some comment6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>Big262Ern\tuclust_test_seqs_7 some comment7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>Big263Ern\tuclust_test_seqs_8 some comment8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>Big264Ern\tuclust_test_seqs_9 some comment9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA""".split('\n')

expected_clusters_w_abundance_default_settings = """>Cluster1;size=326
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT""".split('\n')

expected_clusters_w_abundance_low_setting = """>Cluster1;size=304
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Cluster11;size=22
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCTAACCC
CCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT""".split('\n')

expected_reference_filtered_seqs = """>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTG
ACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGT
>Solemya seq
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAG
TGGCGGACGGGTGAGTA
>usearch_ecoli_seq2
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTG
ACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTCCAT
>Solemya_seq2
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAG
TGGCGGACGGGTGAGTATCAAG""".split('\n')

expected_de_novo_chimeras_default = """>Cluster1;size=52
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACAC
CGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGT
TATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster0;size=50
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCC
CGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCA
TCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGG
>Cluster2;size=45
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCC
CGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>Cluster10;size=43
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCC
CGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTC
CTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>Cluster4;size=40
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>Cluster6;size=40
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTA
CGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGG
>Cluster3;size=30
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCT
TACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTG
TTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>Cluster12;size=19
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTGGTCTTGGTGAGCCGTTACCC
CACCAACTAACTAATACAACGCATGCCCATCCATTACCACCGGAGTTTTCAACCCAAGAAGATGCCTCCCTGGATGTTAT
GGGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTAATGGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTC
GCCGACAAT
>Cluster29;size=18
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGACGCCGCGTC
ACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
TCGCCGG
>Cluster30;size=18
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCC
CGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAAGAGATGCCTCTCGGAGTTAT
TACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTACGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster16;size=16
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCC
CTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTA
TGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACT
>Cluster221;size=1
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTGCCC
CGCCAACTACCTAATCGGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGATTATCTCACCATGCGGCAAAATAATGT
CATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGCGTT
>Cluster218;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGTCACCTTGGTAGGCCGTTACCC
CACCAACTAGCTAATCAGACGCAAGCCCATCTATCAGCGGATTGCTCCTTTTCTAGCTATATCATGCGATACTACTAGCT
TATGCGGTATTAGCAATGATTTCTCACTGTTATTCCCCTCTGATAGGCAGG
>Cluster217;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGATATCCAGAGCCATGCGACCCAGATATATT
ATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCATCCAAGGCAGGTT
>Cluster216;size=1
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGTCAGCTACTGATCGTCGCCTTGGTAGGCCATTACCC
TACCAACTAGCTAATCAGACGCGAGGCCATCTCTCAGCGATAAATCTTTGATATATCTGCCATGCGACAAACATATATTA
TGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTT""".split('\n')

expected_de_novo_chimera_filtered_skew11 = """>Cluster1;size=52
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACAC
CGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGT
TATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster0;size=50
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCC
CGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCA
TCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGG
>Cluster2;size=45
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCC
CGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>Cluster10;size=43
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCC
CGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTC
CTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>Cluster4;size=40
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>Cluster6;size=40
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTA
CGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGG
>Cluster3;size=30
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCT
TACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTG
TTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>Cluster12;size=19
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTGGTCTTGGTGAGCCGTTACCC
CACCAACTAACTAATACAACGCATGCCCATCCATTACCACCGGAGTTTTCAACCCAAGAAGATGCCTCCCTGGATGTTAT
GGGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTAATGGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTC
GCCGACAAT
>Cluster29;size=18
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGACGCCGCGTC
ACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
TCGCCGG
>Cluster30;size=18
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCC
CGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAAGAGATGCCTCTCGGAGTTAT
TACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTACGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster16;size=16
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCC
CTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTA
TGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACT
>Cluster522;size=10
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAT
GCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACT
>Cluster221;size=1
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTGCCC
CGCCAACTACCTAATCGGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGATTATCTCACCATGCGGCAAAATAATGT
CATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGCGTT
>Cluster218;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGTCACCTTGGTAGGCCGTTACCC
CACCAACTAGCTAATCAGACGCAAGCCCATCTATCAGCGGATTGCTCCTTTTCTAGCTATATCATGCGATACTACTAGCT
TATGCGGTATTAGCAATGATTTCTCACTGTTATTCCCCTCTGATAGGCAGG
>Cluster217;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGATATCCAGAGCCATGCGACCCAGATATATT
ATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCATCCAAGGCAGGTT
>Cluster216;size=1
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGTCAGCTACTGATCGTCGCCTTGGTAGGCCATTACCC
TACCAACTAGCTAATCAGACGCGAGGCCATCTCTCAGCGATAAATCTTTGATATATCTGCCATGCGACAAACATATATTA
TGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTT""".split('\n')

expected_cluster_err_seqs = """>Cluster0;size=326
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT""".split('\n')

expected_sorted_by_abundance_no_filter = """>Cluster1;size=52
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACAC
CGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGT
TATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster0;size=50
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCC
CGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCA
TCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGG
>Cluster2;size=45
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCC
CGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>Cluster10;size=43
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCC
CGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTC
CTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>Cluster4;size=40
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>Cluster6;size=40
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTA
CGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGG
>Cluster3;size=30
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCT
TACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTG
TTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>Cluster12;size=19
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTGGTCTTGGTGAGCCGTTACCC
CACCAACTAACTAATACAACGCATGCCCATCCATTACCACCGGAGTTTTCAACCCAAGAAGATGCCTCCCTGGATGTTAT
GGGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTAATGGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTC
GCCGACAAT
>Cluster29;size=18
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGACGCCGCGTC
ACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
TCGCCGG
>Cluster30;size=18
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCC
CGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAAGAGATGCCTCTCGGAGTTAT
TACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTACGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster16;size=16
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCC
CTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTA
TGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACT
>Cluster522;size=10
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAT
GCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACT
>Cluster222;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAT
GCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACT
>Cluster221;size=1
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTGCCC
CGCCAACTACCTAATCGGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGATTATCTCACCATGCGGCAAAATAATGT
CATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGCGTT
>Cluster218;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGTCACCTTGGTAGGCCGTTACCC
CACCAACTAGCTAATCAGACGCAAGCCCATCTATCAGCGGATTGCTCCTTTTCTAGCTATATCATGCGATACTACTAGCT
TATGCGGTATTAGCAATGATTTCTCACTGTTATTCCCCTCTGATAGGCAGG
>Cluster217;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGATATCCAGAGCCATGCGACCCAGATATATT
ATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCATCCAAGGCAGGTT
>Cluster216;size=1
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGTCAGCTACTGATCGTCGCCTTGGTAGGCCATTACCC
TACCAACTAGCTAATCAGACGCGAGGCCATCTCTCAGCGATAAATCTTTGATATATCTGCCATGCGACAAACATATATTA
TGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTT""".split('\n')

expected_abundance_sort_filtered = """>Cluster1;size=52
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACAC
CGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGT
TATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCG
GTCGCCGG
>Cluster0;size=50
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCC
CGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCA
TCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGG
>Cluster2;size=45
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCC
CGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>Cluster10;size=43
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCC
CGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTC
CTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>Cluster4;size=40
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCC
CGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTA
TGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>Cluster6;size=40
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
CGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTA
CGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGG""".split('\n')

expected_derep_seqs = """>seq1;size=2
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAA
>seq2;size=2
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC""".split('\n')

expected_derep_seqs_full_len = """>Cluster0;size=1
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAA
>Cluster1;size=2
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
>Cluster2;size=1
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTT""".split('\n')

expected_len_sorted_seqs = """>chimera
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACCCCTAGGGTGGGAATAACCCGGGGAAACCCGGGCTAATACCGAATAAGACCACAGGAGGCGACTCCAGAGGGTCAAAGGGAGCCTTGGCCTCCCCC
>usearch_ecoli_seq2
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTCCAT
>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGT
>Solemya_seq2
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTATCAAG
>Solemya seq
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTA""".split('\n')

expected_combined_dna_seqs_1_seqs_usearch = """>uclust_test_seqs_0 some comment0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_1 some comment1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_2 some comment2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3 some comment3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>uclust_test_seqs_4 some comment4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>uclust_test_seqs_5 some comment4_again
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>uclust_test_seqs_6 some comment6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>uclust_test_seqs_7 some comment7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_8 some comment8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>uclust_test_seqs_9 some comment9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGT
>Solemya seq
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTA
>usearch_ecoli_seq2
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTCCAT
>Solemya_seq2
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTATCAAG
>chimera
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACCCCTAGGGTGGGAATAACCCGGGGAAACCCGGGCTAATACCGAATAAGACCACAGGAGGCGACTCCAGAGGGTCAAAGGGAGCCTTGGCCTCCCCC""".split('\n')

expected_retained_chimeras_union = """>seq1
ACAGGCC
>seq2
ACAGGCCCCC
>seq3
TTATCCATT
>seq4
ACAGGCCCCC
>seq5
TTATCCATT""".split('\n')

expected_retained_chimeras_intersection = """>seq3
TTATCCATT""".split('\n')

expected_usearch61_denovo_uchime_file = """0.0000\tCluster1;size=52\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\tCluster0;size=50\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\tCluster2;size=45\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\tCluster10;size=43\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\tCluster4;size=40\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\tCluster6;size=40\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\tCluster3;size=30\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0263\tCluster12;size=19\tCluster2;size=45\tCluster1;size=52\tCluster1;size=52\t75.6\t73.3\t76.5\t67.3\t75.6\t20\t21\t26\t6\t1\t3\t*\tN
0.0000\tCluster30;size=18\t*\t*\tCluster6;size=40\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0924\tCluster29;size=18\tCluster6;size=40\tCluster1;size=52\tCluster1;size=52\t92.0\t88.6\t89.0\t86.5\t88.7\t7\t0\t0\t12\t7\t14\t3.3\tN
0.0187\tCluster16;size=16\tCluster2;size=45\tCluster4;size=40\tCluster4;size=40\t94.5\t92.3\t94.1\t90.9\t94.0\t2\t1\t0\t9\t4\t7\t0.5\tN
0.4232\tCluster222;size=1\tCluster4;size=40\tCluster2;size=45\tCluster2;size=45\t100.0\t94.1\t97.3\t91.3\t96.8\t7\t1\t0\t13\t0\t0\t3.2\tY
0.0759\tCluster221;size=1\tCluster16;size=16\tCluster1;size=52\tCluster16;size=16\t74.5\t75.9\t67.3\t66.8\t75.4\t15\t0\t5\t16\t19\t32\t*\tN
0.0107\tCluster218;size=1\tCluster2;size=45\tCluster4;size=40\tCluster4;size=40\t81.7\t80.7\t80.7\t90.6\t78.7\t6\t5\t28\t2\t0\t3\t3.0\tN
0.0086\tCluster217;size=1\tCluster4;size=40\tCluster2;size=45\tCluster4;size=40\t83.1\t83.1\t80.7\t90.8\t82.1\t4\t0\t1\t2\t4\t33\t1.0\tN
0.0000\tCluster216;size=1\t*\t*\tCluster16;size=16\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.4232\tCluster522;size=10\tCluster4;size=40\tCluster2;size=45\tCluster2;size=45\t100.0\t94.1\t97.3\t91.3\t96.8\t7\t1\t0\t13\t0\t0\t3.2\tY""".split('\n')

expected_usearch61_ref_uchime_file = """0.0000\tCluster1;size=52\t*\t*\tseq1\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\tCluster0;size=50\t*\t*\tseq2\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.0000\tCluster2;size=45\t*\t*\tseq3\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN
0.1074\tCluster10;size=43\tmixed_seq\tseq1\tseq1\t70.3\t67.0\t65.1\t54.1\t65.7\t11\t0\t1\t31\t27\t33\t4.6\tN
0.6322\tCluster4;size=40\tmixed_seq\tseq3\tseq3\t96.0\t77.6\t92.5\t73.6\t91.0\t6\t0\t0\t38\t2\t5\t5.0\tY
0.1101\tCluster6;size=40\tseq2\tseq1\tseq1\t82.6\t71.3\t85.2\t69.6\t85.1\t12\t19\t16\t25\t0\t4\t*\tN
0.0258\tCluster3;size=30\tmixed_seq\tseq3\tseq3\t71.6\t66.0\t68.0\t71.1\t66.4\t12\t7\t36\t16\t7\t5\t5.3\tN
0.0263\tCluster12;size=19\tseq3\tseq1\tseq1\t75.6\t73.3\t76.5\t67.3\t75.6\t20\t21\t26\t6\t1\t3\t*\tN
0.0530\tCluster30;size=18\tseq2\tseq1\tseq1\t79.6\t68.3\t85.7\t70.4\t85.9\t8\t24\t16\t25\t0\t6\t*\tN
0.0534\tCluster29;size=18\tseq2\tseq1\tseq1\t80.9\t70.4\t88.3\t70.0\t88.7\t7\t25\t17\t23\t0\t2\t*\tN
0.0699\tCluster16;size=16\tmixed_seq\tseq3\tseq3\t94.0\t74.6\t93.5\t73.6\t91.9\t2\t2\t2\t41\t3\t5\t2.1\tN
1.2277\tCluster222;size=1\tmixed_seq\tseq3\tseq3\t100.0\t78.4\t97.1\t75.5\t96.8\t6\t1\t0\t44\t0\t0\t3.2\tY
0.0855\tCluster221;size=1\tseq3\tseq1\tseq3\t75.8\t77.2\t68.8\t65.1\t72.9\t14\t0\t4\t17\t18\t28\t2.9\tN
0.0174\tCluster218;size=1\tmixed_seq\tseq3\tseq3\t81.7\t70.3\t80.7\t70.3\t78.0\t1\t0\t4\t34\t12\t21\t3.6\tN
0.0713\tCluster217;size=1\tmixed_seq\tseq3\tseq3\t83.3\t77.5\t79.9\t68.6\t79.7\t4\t0\t1\t27\t12\t17\t3.6\tN
0.0505\tCluster216;size=1\tmixed_seq\tseq3\tseq3\t77.5\t72.5\t71.6\t70.1\t72.0\t14\t4\t27\t15\t5\t8\t5.4\tN
1.2277\tCluster522;size=10\tmixed_seq\tseq3\tseq3\t100.0\t78.4\t97.1\t75.5\t96.8\t6\t1\t0\t44\t0\t0\t3.2\tY""".split('\n')

usearch61_dereplicated_uc_lines = """S	0	80	*	*	*	*	*	seq2	*
H	0	80	100.0	*	0	0	*	seq3	seq2
H	0	80	100.0	*	0	0	*	seq4	seq2
S	1	80	*	*	*	*	*	seq1	*
C	0	3	*	*	*	*	*	seq2	*
C	1	1	*	*	*	*	*	seq1	*""".split('\n')

usearch61_clustered_uc_lines = """S	0	80	*	*	*	*	*	seq2;size=3;	*
S	1	80	*	*	*	*	*	seq1;size=1;	*
C	0	1	*	*	*	*	*	seq2;size=3;	*
C	1	1	*	*	*	*	*	seq1;size=1;	*""".split('\n')

usearch61_clustered_uc_lines_ref = """H	3	80	100.0	+	0	0	80M	seq2;size=3;	seq4
H	0	80	100.0	+	0	0	80M	seq1;size=1;	seq1""".split('\n')

usearch61_clustered_ref_lines = """H	0	80	100.0	+	0	0	80M	seq2;size=3;	seq2
N	*	*	*	.	*	*	*	seq1;size=1;	*""".split('\n')

if __name__ == '__main__':
    main()
