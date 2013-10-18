#!/usr/bin/env python

"""
 : provides unit tests for the uclust.py module

Modified from Daniel McDonald's test_cd_hit.py code on Feb-4-2010 """

from os import getcwd, rmdir, remove
from subprocess import Popen, PIPE, STDOUT
from os.path import isfile
from cogent.util.misc import remove_files
from cogent.core.moltype import DNA
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationError, get_tmp_filename

from qiime.pycogent_backports.uclust import (Uclust, 
 uclust_fasta_sort_from_filepath,
 uclust_cluster_from_sorted_fasta_filepath,
 get_output_filepaths,clusters_from_uc_file,
 get_clusters_from_fasta_filepath,
 uclust_search_and_align_from_fasta_filepath,
 process_uclust_pw_alignment_results, UclustParseError)

__author__ = "William Walters"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald","William Walters","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Production"

class UclustTests(TestCase):

    def setUp(self):
        
        self.tmp_unsorted_fasta_filepath = \
         get_tmp_filename(prefix="uclust_test", suffix=".fasta")
        tmp_unsorted_fasta = open(self.tmp_unsorted_fasta_filepath,"w")
        tmp_unsorted_fasta.write('\n'.join(raw_dna_seqs))
        tmp_unsorted_fasta.close()
        
        self.tmp_sorted_fasta_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".fasta")
        tmp_sorted_fasta = open(self.tmp_sorted_fasta_filepath,"w")
        tmp_sorted_fasta.write('\n'.join(sorted_dna_seqs))
        tmp_sorted_fasta.close()
        
        self.tmp_uc_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".uc")
        tmp_uc = open(self.tmp_uc_filepath,"w")
        tmp_uc.write('\n'.join(uc_dna_clusters))
        tmp_uc.close()
         
        self.tmp_clstr_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".clstr")
         
        self.WorkingDir = '/tmp/uclust_test'
        self.tmpdir = '/tmp/'
        
        self.files_to_remove = [self.tmp_unsorted_fasta_filepath,
                                self.tmp_sorted_fasta_filepath,
                                self.tmp_uc_filepath,
                                self.tmp_clstr_filepath]
    
    def tearDown(self):
        remove_files(self.files_to_remove,error_on_missing=False)
        
    def test_fasta_sorting(self):
        """ Should sort fasta seqs from largest to smallest in outfile 
        
        Since a fasta file has to be passed to the app controller for uclust,
        a temporary fasta file is created, and the raw fasta seqs supplied
        in this module are written to it.  This file is sent to the app 
        controller, and the resulting sorted file is compared to the expected
        results to ensure proper function of uclust as called by this app
        controller."""

        test_app = Uclust({'--tmpdir':self.tmpdir})

        
        test_app_res = test_app(data = \
         {'--mergesort':self.tmp_unsorted_fasta_filepath,\
         '--output':self.tmp_sorted_fasta_filepath})

        sorted_fasta_actual = [l.strip() 
            for l in open(test_app_res['Output'].name,"U")]
        sorted_fasta_expected = [l.strip() for l in sorted_dna_seqs if l]
            
        self.assertEqual(sorted_fasta_actual,sorted_fasta_expected)
        
        test_app_res.cleanUp()
    
    def test_parameter_availability(self):
        """ Often used parameters are accessible 
        
            This is just some basic sanity checking.
        
        """
        a = Uclust()
        # if a parameter is not accessible, trying to turn it on will 
        # raise a KeyError
        a.Parameters['--allhits'].on()
        a.Parameters['--libonly'].on()
        a.Parameters['--maxaccepts'].on(42)
        a.Parameters['--maxrejects'].on(42)
        a.Parameters['--rev'].on()

    def test_clustering_fasta_filepath(self):
        """ Should create clusters in uclust format from sorted fasta file 
        
        Since a fasta file has to be passed to the app controller for uclust,
        a temporary fasta file is created, and the sorted seqs supplied
        in this module are written to it.  This file is sent to the app 
        controller, and the resulting uclust file is compared to the expected
        results to ensure proper function of uclust as called by this app
        controller."""
        
        

        test_app = Uclust({'--id':0.9},HALT_EXEC=False)
        test_app_res = test_app(data = \
         {'--input':self.tmp_sorted_fasta_filepath,\
         '--uc':self.tmp_uc_filepath})
        
        uc_file = open(test_app_res['ClusterFile'].name,"U")
        # compare the actual and expect uc files, ignoring comment lines
        uc_file_actual = [l.strip() for l in uc_file 
         if not l.startswith('#')]
        uc_file_expected = [l.strip() for l in uc_dna_clusters 
         if not l.startswith('#')]
        
        self.assertEqual(uc_file_actual, uc_file_expected)
    
        test_app_res.cleanUp()
    
class UclustConvenienceWrappers(TestCase):
    """ Unit tests for uclust convenience wrappers """

    def setUp(self):
        
        self.tmp_unsorted_fasta_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".fasta")
        tmp_unsorted_fasta = open(self.tmp_unsorted_fasta_filepath,"w")
        tmp_unsorted_fasta.write('\n'.join(raw_dna_seqs))
        tmp_unsorted_fasta.close()
        
        self.tmp_raw_dna_seqs_rc_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".fasta")
        tmp_rc_fasta = open(self.tmp_raw_dna_seqs_rc_filepath,"w")
        tmp_rc_fasta.write('\n'.join(raw_dna_seqs_rc))
        tmp_rc_fasta.close()
        
        self.tmp_sorted_fasta_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".fasta")
        tmp_sorted_fasta = open(self.tmp_sorted_fasta_filepath,"w")
        tmp_sorted_fasta.write('\n'.join(sorted_dna_seqs))
        tmp_sorted_fasta.close()
        
        self.tmp_uc_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".uc")
        tmp_uc = open(self.tmp_uc_filepath,"w")
        tmp_uc.write('\n'.join(uc_dna_clusters))
        tmp_uc.close()
        
        self.tmp_clstr_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".clstr")
        
        self.search_align_out1_expected = search_align_out1_expected
        self.search_align_out_fasta_pairs1 = search_align_out_fasta_pairs1
        self.search_align_out_uc1 = search_align_out_uc1
        self.search_align_query1_fp = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".clstr")
        open(self.search_align_query1_fp,'w').write(search_align_query1)
        self.search_align_template1_fp = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".clstr")
        open(self.search_align_template1_fp,'w').write(search_align_template1)
        
        
        self.search_align_out2_expected = search_align_out2_expected
        self.search_align_query2_fp = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".clstr")
        open(self.search_align_query2_fp,'w').write(search_align_query2)
        self.search_align_template2_fp = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".clstr")
        open(self.search_align_template2_fp,'w').write(search_align_template2)
        
        self.ref_dna_seqs_fp = get_tmp_filename(prefix = "uclust_test", suffix = ".fasta")
        open(self.ref_dna_seqs_fp,'w').write(ref_dna_seqs)
        
        self.files_to_remove = [self.tmp_unsorted_fasta_filepath,
                                self.tmp_raw_dna_seqs_rc_filepath,
                                self.tmp_sorted_fasta_filepath,
                                self.tmp_uc_filepath,
                                self.tmp_clstr_filepath,
                                self.search_align_query1_fp,
                                self.search_align_template1_fp,
                                self.search_align_query2_fp,
                                self.search_align_template2_fp,
                                self.ref_dna_seqs_fp]
        
        self.ref_test_clusters1 = ref_test_clusters1
        self.ref_test_failures1 = ref_test_failures1
        self.ref_test_new_seeds1 = ref_test_new_seeds1
        self.ref_test_clusters2 = ref_test_clusters2
        self.ref_test_failures2 = ref_test_failures2
        self.ref_test_new_seeds2 = ref_test_new_seeds2
        self.uc_dna_clusters = uc_dna_clusters
        self.uc_lines1 = uc_lines1
        self.uc_lines_w_multiple_hits_per_query = \
         uc_lines_w_multiple_hits_per_query
        self.uc_lines_overlapping_lib_input_seq_ids = \
         uc_lines_overlapping_lib_input_seq_ids
        
    def tearDown(self):
        remove_files(self.files_to_remove,error_on_missing=False)


    def test_uclust_fasta_sort_from_filepath(self):
        """ Given an unsorted fasta filepath, will return sorted file """
        
        app_res = \
         uclust_fasta_sort_from_filepath(self.tmp_unsorted_fasta_filepath)
        
        sorted_fasta_actual = [l.strip() 
            for l in open(app_res['Output'].name,"U")]
        sorted_fasta_expected = [l.strip() for l in sorted_dna_seqs if l]
        
        self.assertEqual(sorted_fasta_actual,sorted_fasta_expected)
        
        app_res.cleanUp()
        
        
    def test_clusters_from_uc_file(self):
        """ clusters_from_uc_file functions as expected """

        expected_clusters = {'s2':['s2','s3']}
        expected_failures = ['s1']
        expected_new_seeds = ['s2']
        self.assertEqual(clusters_from_uc_file(self.uc_lines1),
         (expected_clusters,expected_failures,expected_new_seeds))
    
    def test_clusters_from_uc_file_multiple_hits(self):
        """ clusters_from_uc_file handles error_on_multiple_hits correctly
        """
        # when a query hits multiple hits and error_on_multiple_hits=True
        # an error should be raised
        self.assertRaises(UclustParseError,
                          clusters_from_uc_file,
                          self.uc_lines_w_multiple_hits_per_query,
                          error_on_multiple_hits=True)
        
        # when a query hits multiple hits and error_on_multiple_hits=False
        # the query should show up in multiple clusters
        actual = clusters_from_uc_file(self.uc_lines_w_multiple_hits_per_query,
                                       error_on_multiple_hits=False)
        expected_clusters = {'s2':['s2','s3'],
                             's4':['s4','s3']}
        expected_failures = ['s1']
        expected_new_seeds = ['s2','s4']
        self.assertEqual(actual,
         (expected_clusters,expected_failures,expected_new_seeds))

    def test_clusters_from_uc_file_error(self):
        """ clusters_from_uc_file raises error when lib/input seq ids overlap
        """
        self.assertRaises(UclustParseError,
                          clusters_from_uc_file,
                          self.uc_lines_overlapping_lib_input_seq_ids)
        
        
    def test_uclust_cluster_from_sorted_fasta_filepath(self):
        """ Given a sorted fasta filepath, will return uclust (.uc) file """
        

        app_res = \
         uclust_cluster_from_sorted_fasta_filepath(self.tmp_sorted_fasta_filepath, \
         percent_ID = 0.90,HALT_EXEC=False)

        
        uc_file = open(app_res['ClusterFile'].name,"U")
        # compare the actual and expect uc files, ignoring comment lines
        uc_file_actual = [l.strip() for l in uc_file 
         if not l.startswith('#')]
        uc_file_expected = [l.strip() for l in uc_dna_clusters 
         if not l.startswith('#')]
        
        self.assertEqual(uc_file_actual, uc_file_expected)
        app_res.cleanUp()
    
    def test_get_output_filepaths(self):
        """ Properly generates output filepath names """
        
        uc_res = \
         get_output_filepaths("/tmp/","test_seqs.fasta")
        
        self.assertEqual(uc_res, "/tmp/test_seqs_clusters.uc")

        

        
    def test_get_clusters_from_fasta_filepath(self):
        """ Tests for return of lists of OTUs from given fasta filepath """
        
        clusters_res = \
         get_clusters_from_fasta_filepath(self.tmp_unsorted_fasta_filepath, \
          original_fasta_path = None, percent_ID = 0.90, save_uc_files=False)
        expected_cluster_list.sort()
        expected_failure_list.sort()
        expected_new_seed_list.sort()
        clusters_res[0].sort()
        clusters_res[1].sort()
        clusters_res[2].sort()
        self.assertEqual(clusters_res,(expected_cluster_list,
                                       expected_failure_list,
                                       expected_new_seed_list))
                                       
    def test_get_clusters_from_fasta_filepath_reference_db_only(self):
        """ Correct clusters returned when clustering against a database only 
        """
        clusters_res = get_clusters_from_fasta_filepath(
          self.tmp_unsorted_fasta_filepath,
          original_fasta_path = None, 
          save_uc_files=False,
          max_accepts=7,max_rejects=12,
          percent_ID = 0.90,
          subject_fasta_filepath=self.ref_dna_seqs_fp,
          suppress_new_clusters=True,
          HALT_EXEC=False)
        
        self.ref_test_clusters1.sort()
        self.ref_test_failures1.sort()
        self.ref_test_new_seeds1.sort()
        
        clusters_res[0].sort()
        clusters_res[1].sort()
        clusters_res[2].sort()
        self.assertEqual(clusters_res,(self.ref_test_clusters1,
                                       self.ref_test_failures1,
                                       self.ref_test_new_seeds1))
                                       
    def test_get_clusters_from_fasta_filepath_extending_reference_db(self):
        """ Correct clusters when clustering against db and adding new clusters
        """
        clusters_res = get_clusters_from_fasta_filepath(
          self.tmp_unsorted_fasta_filepath,
          original_fasta_path = None,
          max_accepts=7,max_rejects=12,
          percent_ID = 0.90,
          subject_fasta_filepath=self.ref_dna_seqs_fp,
          suppress_new_clusters=False,enable_rev_strand_matching=True,
          HALT_EXEC=False,
          save_uc_files=False)
        
        self.ref_test_clusters2.sort()
        self.ref_test_failures2.sort()
        self.ref_test_new_seeds2.sort()
        
        clusters_res[0].sort()
        clusters_res[1].sort()
        clusters_res[2].sort()
        self.assertEqual(clusters_res,(self.ref_test_clusters2,
                                       self.ref_test_failures2,
                                       self.ref_test_new_seeds2))
        

    def test_get_clusters_from_fasta_filepath_optimal(self):
        """ Test OTUs from filepath functions with optimal
        """
        # need to compile a small test where optimal has an affect --
        # this currently is only testing that we don't get a failure with
        # optimal
        clusters_res = \
         get_clusters_from_fasta_filepath(self.tmp_unsorted_fasta_filepath,
          original_fasta_path = None, save_uc_files=False,
          percent_ID = 0.90, optimal = True)
        expected_cluster_list.sort()
        expected_failure_list.sort()
        expected_new_seed_list.sort()
        clusters_res[0].sort()
        clusters_res[1].sort()
        clusters_res[2].sort()
        
        self.assertEqual(clusters_res,(expected_cluster_list,
                                       expected_failure_list,
                                       expected_new_seed_list))

        
    def test_get_clusters_from_fasta_filepath_suppress_sort(self):
        """ Test OTUs from filepath functions with suppress sort
        """
        expected = [['uclust_test_seqs_0'], ['uclust_test_seqs_1'],
                    ['uclust_test_seqs_2'], ['uclust_test_seqs_3'],
                    ['uclust_test_seqs_4'], ['uclust_test_seqs_5'],
                    ['uclust_test_seqs_6', 'uclust_test_seqs_8'],
                    ['uclust_test_seqs_7'], ['uclust_test_seqs_9']]
        clusters_res = \
         get_clusters_from_fasta_filepath(self.tmp_unsorted_fasta_filepath,
          original_fasta_path = None,
          percent_ID = 0.90, suppress_sort = True, save_uc_files=False)
        expected_cluster_list.sort()
        expected_failure_list.sort()
        expected_new_seed_list.sort()
        clusters_res[0].sort()
        clusters_res[1].sort()
        clusters_res[2].sort()
        
        self.assertEqual(clusters_res,(expected_cluster_list,
                                       expected_failure_list,
                                       expected_new_seed_list))
        
    def test_get_clusters_from_fasta_filepath_rev_strand_match(self):
        """ Test OTUs from filepath functions with rev strand match
        """
        # seq and its rc don't cluster when enable_rev_strand_matching = False
        expected_cluster_list = [['uclust_test_seqs_0'], ['uclust_test_seqs_0_rc']]
        expected_failure_list = []
        expected_new_seed_list = ['uclust_test_seqs_0', 'uclust_test_seqs_0_rc']
        clusters_res = \
         get_clusters_from_fasta_filepath(self.tmp_raw_dna_seqs_rc_filepath,
          original_fasta_path = None, save_uc_files=False,
          percent_ID = 0.90, enable_rev_strand_matching = False)
        
        expected_cluster_list.sort()
        expected_failure_list.sort()
        expected_new_seed_list.sort()
        clusters_res[0].sort()
        clusters_res[1].sort()
        clusters_res[2].sort()
        self.assertEqual(clusters_res,(expected_cluster_list,
                                       expected_failure_list,
                                       expected_new_seed_list))
        
        # seq and its rc cluster when enable_rev_strand_matching = False
        expected_cluster_list = [['uclust_test_seqs_0', 'uclust_test_seqs_0_rc']]
        expected_failure_list = []
        expected_new_seed_list = ['uclust_test_seqs_0']
        clusters_res = \
         get_clusters_from_fasta_filepath(self.tmp_raw_dna_seqs_rc_filepath,
          original_fasta_path = None, save_uc_files=False,
          percent_ID = 0.90, enable_rev_strand_matching = True)
        
        expected_cluster_list.sort()
        expected_failure_list.sort()
        expected_new_seed_list.sort()
        clusters_res[0].sort()
        clusters_res[1].sort()
        clusters_res[2].sort()
        self.assertEqual(clusters_res,(expected_cluster_list,
                                       expected_failure_list,
                                       expected_new_seed_list))
        
    def test_process_uclust_pw_alignment_results(self):
        """parsing of pairwise alignment fasta pairs file functions as expected
        """
        actual = list(process_uclust_pw_alignment_results(\
         self.search_align_out_fasta_pairs1,self.search_align_out_uc1))
        expected = self.search_align_out1_expected
        
        # iterate over results so error output will highlight the bad match
        for a,e in zip(actual,expected):
            self.assertEqual(a,e)
        
        # make sure the full result objects are the same
        self.assertEqual(actual,expected)
        
    def test_uclust_search_and_align_from_fasta_filepath(self):
        """ uclust_search_and_align_from_fasta_filepath functions as expected """
        # rev comp matches allowed (default)
        actual = list(uclust_search_and_align_from_fasta_filepath(
         self.search_align_query1_fp,self.search_align_template1_fp))
        self.assertEqual(actual,self.search_align_out1_expected)
        
        # rev comp matches not allowed
        actual = list(uclust_search_and_align_from_fasta_filepath(
         self.search_align_query1_fp,self.search_align_template1_fp,
         enable_rev_strand_matching=False))
        self.assertEqual(actual,self.search_align_out1_expected[:2])
        
    def test_uclust_search_and_align_from_fasta_filepath_protein(self):
        """ uclust_search_and_align_from_fasta_filepath functions with protein """
        # rev comp matches allowed (default)
        actual = list(uclust_search_and_align_from_fasta_filepath(
         self.search_align_query2_fp,self.search_align_template2_fp))
        self.assertEqual(actual,self.search_align_out2_expected)
        
    def test_uclust_supported_version(self):
        """uclust version is supported """
        command = 'uclust --version'
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split('v')[-1].strip('q')
        try:
            version = tuple(map(int,version_string.split('.')))
            acceptable_version = version >= (1,2,22)
        except ValueError:
            acceptable_version = False
        
        self.assertTrue(acceptable_version,\
         "Unsupported uclust version. 1.2.22 or later "+\
         "is required, but running %s." % version_string)

raw_dna_seqs = """>uclust_test_seqs_0
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_1
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>uclust_test_seqs_2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>uclust_test_seqs_4
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_5
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>uclust_test_seqs_6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>uclust_test_seqs_7
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>uclust_test_seqs_9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
""".split('\n')

ref_dna_seqs = """>ref1 25 random bases appended to uclust_test_seqs_0 and one mismatch
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATATTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCTATAGCAGCCCCAGCGTTTACTTCTA
>ref2 15 random bases prepended to uclust_test_seqs_1 and one mismatch
GCTGCGGCGTCCTGCGCCACGGTGGGTACAACACGTCCACTACATCTGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>ref3 5 random bases prepended and 10 random bases appended to uclust_test_seqs_2
ATAGGCCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACTGCCTGATTCA
>ref4 exact match to uclust_test_seqs_3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
"""

ref_test_clusters1 = [['uclust_test_seqs_0'],['uclust_test_seqs_1'],
                      ['uclust_test_seqs_2'],['uclust_test_seqs_3']]
ref_test_failures1 = ['uclust_test_seqs_4','uclust_test_seqs_5',
                      'uclust_test_seqs_6','uclust_test_seqs_7',
                      'uclust_test_seqs_8','uclust_test_seqs_9']
ref_test_new_seeds1 = []

ref_test_clusters2 = [['uclust_test_seqs_0'],['uclust_test_seqs_1'],
                      ['uclust_test_seqs_2'],['uclust_test_seqs_3'],
                      ['uclust_test_seqs_4'],['uclust_test_seqs_5'],
                      ['uclust_test_seqs_6','uclust_test_seqs_8'],
                      ['uclust_test_seqs_7'],['uclust_test_seqs_9']]
ref_test_failures2 = []
ref_test_new_seeds2 = ['uclust_test_seqs_4','uclust_test_seqs_5','uclust_test_seqs_6',
 'uclust_test_seqs_7','uclust_test_seqs_9']


raw_dna_seqs_rc = """>uclust_test_seqs_0
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_0_rc
AGCTCTGACACAAAACTGACGTGATGTGCCTTAAGTATCCAACCCGTTGGATGGGACGTCTTGTAGCCACCGT
""".split('\n')

sorted_dna_seqs=""">uclust_test_seqs_7
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_4
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>uclust_test_seqs_1
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>uclust_test_seqs_5
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>uclust_test_seqs_6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>uclust_test_seqs_0
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>uclust_test_seqs_9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
""".split('\n')

# Clusters are created at a 0.90% identity
uc_dna_clusters= """# uclust --input /tmp/uclust_testBGwZvcikrbNefYGRTk0u.fasta --id 0.9 --uc /tmp/uclust_testrbcO0CyBVpV9AwH3OIK1.uc
# version=1.1.577
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	80	*	*	*	*	*	uclust_test_seqs_7	*
S	1	79	*	*	*	*	*	uclust_test_seqs_4	*
S	2	78	*	*	*	*	*	uclust_test_seqs_2	*
S	3	77	*	*	*	*	*	uclust_test_seqs_3	*
S	4	76	*	*	*	*	*	uclust_test_seqs_1	*
S	5	75	*	*	*	*	*	uclust_test_seqs_5	*
S	6	74	*	*	*	*	*	uclust_test_seqs_6	*
S	7	73	*	*	*	*	*	uclust_test_seqs_0	*
H	6	72	91.7	+	0	0	2I72M	uclust_test_seqs_8	uclust_test_seqs_6
S	8	71	*	*	*	*	*	uclust_test_seqs_9	*
C	0	1	*	*	*	*	*	uclust_test_seqs_7	*
C	1	1	*	*	*	*	*	uclust_test_seqs_4	*
C	2	1	*	*	*	*	*	uclust_test_seqs_2	*
C	3	1	*	*	*	*	*	uclust_test_seqs_3	*
C	4	1	*	*	*	*	*	uclust_test_seqs_1	*
C	5	1	*	*	*	*	*	uclust_test_seqs_5	*
C	6	2	91.7	*	*	*	*	uclust_test_seqs_6	*
C	7	1	*	*	*	*	*	uclust_test_seqs_0	*
C	8	1	*	*	*	*	*	uclust_test_seqs_9	*""".split('\n')

expected_cluster_list=[['uclust_test_seqs_7'], ['uclust_test_seqs_4'], ['uclust_test_seqs_2'], ['uclust_test_seqs_3'], ['uclust_test_seqs_1'], ['uclust_test_seqs_5'], ['uclust_test_seqs_6', 'uclust_test_seqs_8'], ['uclust_test_seqs_0'], ['uclust_test_seqs_9']]
expected_failure_list = []
expected_new_seed_list = ['uclust_test_seqs_7', 'uclust_test_seqs_4', 'uclust_test_seqs_2',
 'uclust_test_seqs_3', 'uclust_test_seqs_1', 'uclust_test_seqs_5', 'uclust_test_seqs_6',
 'uclust_test_seqs_0', 'uclust_test_seqs_9']

search_align_query1 = """>1_like
TACGGCTACCTTGTTACGACTTCATCCCAATCATTTGTTCCACCTTCGACGGCTA
>2_like
ATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAG
>2_like_rc
CTTAGTTGCCATCCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCAT
>rand
TTGCGACGAGCGGACGGCCGGGTGTATGTCGTCATATATATGTGTCTGCCTATCGTTACGTACACTCGTCGTCT
"""

search_align_template1 = """>1
AGAAAGGAGGTGATCCAGCCGCACCTTCCGATACGGCTACCTTGTTACGACTTCACCCCAATCATTTGTTCCACCTTCGACGGCTAGCTCCAAATGGTTACTCCACCGGCTTCGGGTGTTACAAACTC
>2
AGCCCAAATCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAGCTTAAGGGTTGCGCT
"""

search_align_query2 = """>1_like
PRTEINACYYPL
>2_like
AGGYTPPLVN
>rand
GGTYPARREE
"""

search_align_template2 = """>1
PRTELNACYYPL
>2
AGGYTRPPLVN
"""

search_align_out2_expected = [
 ('1_like','1','PRTEINACYYPL','PRTELNACYYPL',91.70000),
 ('2_like','2','AGGYT-PPLVN','AGGYTRPPLVN',100.0)]

search_align_out_fasta_pairs1 = """>1_like
-------------------------------TACGGCTACCTTGTTACGACTTCATCCCAATCATTTGTTCCACCTTCGACGGCTA------------------------------------------
>1+
AGAAAGGAGGTGATCCAGCCGCACCTTCCGATACGGCTACCTTGTTACGACTTCACCCCAATCATTTGTTCCACCTTCGACGGCTAGCTCCAAATGGTTACTCCACCGGCTTCGGGTGTTACAAACTC

>2_like
-------------------ATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAG---------------
>2+
AGCCCAAATCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAGCTTAAGGGTTGCGCT

>2_like_rc
---------------CTTAGTTGCCATCCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCAT-------------------
>2-
AGCGCAACCCTTAAGCTTAGTTGCCATCCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGATTTGGGCT
""".split('\n')

search_align_out_uc1 = """# uclust --input sm_query.fasta --lib sm_template.fasta --id 0.75 --libonly --rev --maxaccepts 8 --maxrejects 32 --fastapairs sm_pw.fasta --uc sm_result.uc
# version=1.1.577
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
L	0	128	*	*	*	*	*	1	*
H	0	55	98.2	+	0	0	31I55M42I	1_like	1
L	1	92	*	*	*	*	*	2	*
H	1	58	100.0	+	0	0	19I58M15I	2_like	2
H	1	58	100.0	-	0	0	15I58M19I	2_like_rc	2
N	*	74	*	*	*	*	*	rand	*
D	0	2	*	*	*	*	98.2	1	*
D	1	3	*	*	*	*	100.0	2	*
""".split('\n')

search_align_out1_expected = [
         ('1_like','1','-------------------------------TACGGCTACCTTGTTACGACTTCATCCCAATCATTTGTTCCACCTTCGACGGCTA------------------------------------------','AGAAAGGAGGTGATCCAGCCGCACCTTCCGATACGGCTACCTTGTTACGACTTCACCCCAATCATTTGTTCCACCTTCGACGGCTAGCTCCAAATGGTTACTCCACCGGCTTCGGGTGTTACAAACTC',98.2),
         
         ('2_like','2','-------------------ATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAG---------------','AGCCCAAATCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAGCTTAAGGGTTGCGCT',100.0),\
         
         ('2_like_rc RC','2','-------------------ATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAG---------------','AGCCCAAATCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAGCTTAAGGGTTGCGCT',100.0)]

uc_lines1 = """# uclust --input q.fasta --lib r.fasta --uc results.uc --id 0.90 --libonly --rev
# version=1.1.579
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
N	*	80	*	*	*	*	*	s1 some comment	*
S	4	80	*	*	*	*	*	s2 some other comment	*
H	2	78	100.0	+	0	0	5I78M10I	s3 yet another comment	s2""".split('\n')

uc_lines_w_multiple_hits_per_query = """# uclust --input q.fasta --lib r.fasta --uc results.uc --id 0.90 --libonly --rev
# version=1.1.579
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
N	*	80	*	*	*	*	*	s1 some comment	*
S	4	80	*	*	*	*	*	s2 some other comment	*
S	4	80	*	*	*	*	*	s4	*
H	2	78	100.0	+	0	0	5I78M10I	s3 yet another comment	s2
H	2	78	98.0	+	0	0	5I78M10I	s3 yet another comment	s4
""".split('\n')

uc_lines_overlapping_lib_input_seq_ids = """# uclust --maxrejects 32 --input /tmp/OtuPickerbb092OWRWLWqlBR2BmTZ.fasta --id 0.97 --uc /tmp/uclust_clustersLf5Oqv0SvGTZo1mVWBqK.uc --rev --usersort --maxaccepts 8 --lib r.fasta
# version=1.1.16
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	1	24	*	*	*	*	*	3	*
H	1	24	100.0	+	0	0	24M	4	3
L	0	54	*	*	*	*	*	3	*
H	0	54	100.0	+	0	0	54M	2	3
D	0	2	*	*	*	*	100.0	3	*
C	1	2	100.0	*	*	*	*	3	*
""".split('\n')

if __name__ == '__main__':
    main()
