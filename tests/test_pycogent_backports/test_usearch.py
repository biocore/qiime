#!/usr/bin/env python

"""
provides unit tests for the usearch.py module
"""

__author__ = "William Walters"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["William Walters", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Development"

from os.path import isfile, basename
from shutil import rmtree
from glob import glob

from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationError, get_tmp_filename

from qiime.util import create_dir

from qiime.pycogent_backports.usearch import (Usearch,
 clusters_from_blast_uc_file, usearch_fasta_sort_from_filepath,
 usearch_dereplicate_exact_subseqs, usearch_dereplicate_exact_seqs,
 usearch_sort_by_abundance,
 usearch_cluster_error_correction, usearch_chimera_filter_de_novo,
 usearch_chimera_filter_ref_based, usearch_cluster_seqs,
 enumerate_otus, assign_reads_to_otus, usearch_qf, concatenate_fastas,
 get_retained_chimeras)

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
          self.tmp_dna_seqs_ref_otu_picking]
          
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
        expected_clusters =  {'1': ['Solemya', 'Solemya_seq2'], 
        '0': ['usearch_ecoli_seq', 'usearch_ecoli_seq2']}
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

dna_seqs_with_dups=""">seq1
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTAA
>seq2
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
>seq3
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCC
>seq4
GCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTT"""

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

expected_cluster_err_seqs = """>Cluster0;size=2
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


if __name__ == '__main__':
    main()
