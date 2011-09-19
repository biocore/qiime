#!/usr/bin/env python

"""
provides unit tests for the usearch.py module
"""


from os.path import isfile, basename
from shutil import rmtree
from glob import glob

from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationError, get_tmp_filename

from qiime.util import create_dir

from qiime.pycogent_backports.usearch import (Usearch,
 clusters_from_blast_uc_file, usearch_fasta_sort_from_filepath,
 usearch_dereplicate_exact_subseqs, usearch_sort_by_abundance,
 usearch_cluster_error_correction, usearch_chimera_filter_de_novo,
 usearch_chimera_filter_ref_based, usearch_cluster_seqs,
 enumerate_otus, assign_reads_to_otus, otu_pipe)

__author__ = "William Walters"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Development"

class UsearchTests(TestCase):

    def setUp(self):
        # create the temporary input files
        self.dna_seqs_1 = dna_seqs_1
        self.dna_seqs_2 = dna_seqs_usearch
        self.ref_database = usearch_ref_seqs1
        self.dna_seqs_with_abundance = dna_seqs_with_abundance
        self.de_novo_chimera_seqs = de_novo_chimera_seqs
        self.dna_seqs_with_dups = dna_seqs_with_dups
        
        
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
        
       
        self._files_to_remove =\
         [self.tmp_seq_filepath1, self.tmp_seq_filepath2,
          self.tmp_ref_database, self.tmp_seqs_w_abundance,
          self.tmp_de_novo_chimera_seqs, self.tmp_dna_seqs_with_dups]
          
        self._dirs_to_remove = []
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        if self._dirs_to_remove:
            for curr_dir in self._dirs_to_remove:
                rmtree(curr_dir)
            
    def test_otu_pipe(self):
        """ Main program loop test, with default parameters """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = otu_pipe(self.tmp_seq_filepath2,
                                      output_dir = self.tmp_dir,
                                      db_filepath = self.tmp_ref_database,
                                      minsize = 1,
                                      remove_usearch_logs=True)
                                      
        expected_clusters = {'1': ['Solemya', 'Solemya_seq2'],
         '0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2']}
        expected_failures = ['chimera']
        
        self.assertEqual(clusters, expected_clusters)
        self.assertEqual(failures, expected_failures)
        
        
    def test_otu_pipe_no_ref_database(self):
        """ Main program loop with no reference chimera testing """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = otu_pipe(self.tmp_seq_filepath2,
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

        
    def test_otu_pipe_disabled_filters(self):
        """ Returns expected clustering with no filtering """
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = otu_pipe(self.tmp_seq_filepath2,
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
        
    def test_otu_pipe_generates_logs(self):
        """ Generates expected log files """
        
        curr_output_dir = get_tmp_filename(prefix='/UsearchLogTest_',suffix='/')
         
        create_dir(curr_output_dir)
        
         
        self._dirs_to_remove.append(curr_output_dir)
        
        # cluster size filtering set to 1 instead of default 4
        clusters, failures = otu_pipe(self.tmp_seq_filepath2,
                                      output_dir = curr_output_dir,
                                      db_filepath = self.tmp_ref_database,
                                      minsize = 1,
                                      remove_usearch_logs=False)
                                      
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

    def test_assign_reads_to_otus(self):
        """ Properly assigns reads back to original ID """
        
        
        app_result, output_filepath =\
         assign_reads_to_otus(original_fasta = self.tmp_ref_database,
                              filtered_fasta = self.tmp_seq_filepath2,
                              remove_usearch_logs = True)
                              
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
                              percent_id = 0.97)
        
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
                              percent_id = 0.999)
        
        self._files_to_remove.append(output_filepath)
                              
        actual_clusters = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_clusters,
         self.expected_clusters_w_abundance_low_setting)
        
    def test_usearch_chimera_filter_ref_based(self):
        """ Properly detects chimeras against reference database """
        
        app_result, output_filepath =\
         usearch_chimera_filter_ref_based(self.tmp_seq_filepath2,
                                          self.tmp_ref_database,
                                          remove_usearch_logs=True)
                                          
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
                                        abundance_skew = 2)
                                        
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = \
         [line.strip() for line in open(output_filepath, "U")]
         
        self.assertEqual(actual_seqs, self.expected_de_novo_chimeras_default)
        
    def test_usearch_chimera_filter_de_novo_abundance_skew(self):
        """ Properly detects de novo chimeras with skew changes """
        
        app_result, output_filepath =\
         usearch_chimera_filter_de_novo(self.tmp_de_novo_chimera_seqs,
                                        remove_usearch_logs=True,
                                        abundance_skew = 11)
                                        
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
                                          percent_id_err = 0.97)
        
        self._files_to_remove.append(output_filepath)
                              
        actual_clusters = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_clusters,
         self.expected_cluster_err_seqs)
        
    def test_usearch_sort_by_abundance(self):
        """ Properly sorts fasta by abundance """
        
        app_result, output_filepath =\
         usearch_sort_by_abundance(self.tmp_de_novo_chimera_seqs,
                                   remove_usearch_logs=True)
                                   
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_seqs,
         self.expected_sorted_by_abundance_no_filter)
         
    def test_usearch_sort_by_abundance_filter(self):
        """ Properly sorts fasta by abundance, filters low count otus """
        
        app_result, output_filepath =\
         usearch_sort_by_abundance(self.tmp_de_novo_chimera_seqs,
                                   remove_usearch_logs=True,
                                   minsize = 40)
                                   
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_seqs,
         self.expected_abundance_sort_filtered)

    def test_usearch_dereplicate_exact_subseqs(self):
        """ Properly dereplicates fasta file """
        
        app_result, output_filepath =\
         usearch_dereplicate_exact_subseqs(self.tmp_dna_seqs_with_dups,
                                           remove_usearch_logs=True)
                                           
        self._files_to_remove.append(output_filepath)
        
        actual_seqs = [line.strip() for line in open(output_filepath, "U")]
        
        self.assertEqual(actual_seqs, self.expected_derep_seqs)
        
    def test_usearch_fasta_sort_from_filepath(self):
        """ Properly sorts fasta according to seq length """
        
        app_result, output_filepath =\
         usearch_fasta_sort_from_filepath(self.tmp_seq_filepath2,
                                          remove_usearch_logs=True)
                                          
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

# Long strings for test files, output, etc.
# *************************************************

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

dna_seqs_with_dups=""">seq1
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAT
>seq2
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
>seq3
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAT"""

# Expected output file data
uc_lines1 = """# usearch --id 0.97 --uc usearch_picked_otus/assign_reads_to_otus.uc --query seqs.fna --global --db usearch_picked_otus/enumerated_otus.fasta
# version=4.2.66
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
H	42	217	99.1	+	0	0	217MI	PC.test2_1 FLP3FBN02xELBSXx orig_bc=ACAGAGTCGGCG new_bc=ACAGAGTCGGCG,FLP3FBN02x bc_diffs=0	42
H	42	217	99.1	+	0	0	217MI	PC.test1_2 FLP3FBN03ELBSXx orig_bc=ACAGAGTCGGCG new_bc=ACAGAGTCGGCG,FLP3FBN03 bc_diffs=0	42
H	42	217	99.1	+	0	0	217MI	PC.634_3 FLP3FBN01ELBSX orig_bc=TCAGAGTCGGCT new_bc=ACAGAGTCGGCT,FLP3FBN01 bc_diffs=1	42
H	19	243	100.0	+	0	0	25MI218M	PC.634_4 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT,FLP3FBN01 bc_diffs=0	19
N	*	219	*	*	*	*	*	PC.481_6 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG,FLP3FBN01 bc_diffs=0	*
H	6	211	99.5	+	0	0	211M	PC.269_5 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA,FLP3FBN01 bc_diffs=0	6
""".split('\n')

expected_otu_assignments = """# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
H	2	199	97.5	+	0	0	119M80D	ref1 ecoli sequence	usearch_ecoli_seq2
N	*	178	*	*	*	*	*	EU199232 1 1236 Bacteria/Deltaproteobacteria/Desulfurella - Hippea/uncultured	*
H	1	180	100.0	+	0	0	97M83D	L07864 1 1200 Bacteria/Beta Gammaproteobacteria/Solemya symbiont	Solemya seq""".split('\n')

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

expected_enumerated_fasta_added_options = """>Big255Ern	uclust_test_seqs_0 some comment0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>Big256Ern	uclust_test_seqs_1 some comment1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>Big257Ern	uclust_test_seqs_2 some comment2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>Big258Ern	uclust_test_seqs_3 some comment3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>Big259Ern	uclust_test_seqs_4 some comment4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>Big260Ern	uclust_test_seqs_5 some comment4_again
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>Big261Ern	uclust_test_seqs_6 some comment6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>Big262Ern	uclust_test_seqs_7 some comment7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>Big263Ern	uclust_test_seqs_8 some comment8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>Big264Ern	uclust_test_seqs_9 some comment9
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
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAT
>seq2;size=1
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC""".split('\n')

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

if __name__ == '__main__':
    main()
