#!/usr/bin/env python

"""
test_uclust.py : provides unit tests for the uclust.py module

Modified from Daniel McDonald's test_cd_hit.py code on Feb-4-2010 """

from os import getcwd, rmdir, remove
from os.path import isfile
from cogent.util.misc import remove_files
from cogent.core.moltype import DNA
from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename, ApplicationError
from pynast.pycogent_backports.uclust import (UclustFastaSort, 
 uclust_fasta_sort_from_filepath,
 UclustCreateClusterFile, uclust_cluster_from_sorted_fasta_filepath,
 UclustConvertToCdhit, uclust_convert_uc_to_cdhit_from_filepath,
 parse_uclust_clstr_file, get_output_filepaths,
 get_clusters_from_fasta_filepath, process_uclust_blast_result,
 uclust_search_and_align_from_fasta_filepath)

__author__ = "William Walters"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Daniel McDonald","William Walters","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0.dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Development"

class UclustFastaSort_Tests(TestCase):
    """ Tests for UclustFastaSort application controller"""
    
    def setUp(self):
        
        self.tmp_unsorted_fasta_filepath = \
         get_tmp_filename(prefix="uclust_test", suffix="fasta")
        self.tmp_sorted_fasta_filepath = get_tmp_filename(prefix="uclust_test",\
         suffix="fasta")
        self.WorkingDir = '/tmp/uclust_test'
        self.tmpdir = '/tmp/'
        
    def tearDown(self):
        if isfile(self.tmp_unsorted_fasta_filepath):
            remove(self.tmp_unsorted_fasta_filepath)
        if isfile(self.tmp_sorted_fasta_filepath):
            remove(self.tmp_sorted_fasta_filepath)

    
    def test_base_command(self):
        """ UclustFastaSort should return the correct BaseCommand """
        c = UclustFastaSort()
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "',getcwd(),'/"; ','uclust']))
        c.Parameters['--mergesort'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "',getcwd(),'/"; ','uclust --mergesort "seq.txt"']))
        c.Parameters['--output'].on('sorted_output.fasta')
        c.Parameters['--tmpdir'].on(self.tmpdir)
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "',getcwd(),'/"; ','uclust --mergesort "seq.txt" '+\
         '--output "sorted_output.fasta" --tmpdir "/tmp/"']))

    def test_changing_working_dir(self):
        """ UclustFastaSort BaseCommand should change according to WorkingDir"""
        
        c = UclustFastaSort(WorkingDir=self.WorkingDir)
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "','/tmp/uclust_test','/"; ','uclust']))
        c = UclustFastaSort()
        c.WorkingDir = '/tmp/uclust_test2'
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "','/tmp/uclust_test2','/"; ','uclust']))
         
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/uclust_test')
        rmdir('/tmp/uclust_test2')
        
    def test_sort_fasta_from_fasta_filepath(self):
        """ Should sort fasta seqs from largest to smallest in outfile 
        
        Since a fasta file has to be passed to the app controller for uclust,
        a temporary fasta file is created, and the raw fasta seqs supplied
        in this module are written to it.  This file is sent to the app 
        controller, and the resulting sorted file is compared to the expected
        results to ensure proper function of uclust as called by this app
        controller."""
        

        tmp_unsorted_fasta = open(self.tmp_unsorted_fasta_filepath,"w")
        for line in raw_dna_seqs:
            tmp_unsorted_fasta.write(line)

        tmp_unsorted_fasta.close()

        test_app = UclustFastaSort()

        
        test_app_res = test_app(data = \
         {'--mergesort':self.tmp_unsorted_fasta_filepath,\
         '--output':self.tmp_sorted_fasta_filepath,\
         '--tmpdir':self.tmpdir})

        sorted_fasta = open(test_app_res['SortedFasta'].name,"U")
        sorted_fasta_res = []
        for line in sorted_fasta:
            sorted_fasta_res.append(line)
            
        self.assertEqual(sorted_fasta_res, sorted_dna_seqs)
        
        test_app_res.cleanUp()

class UclustCreateClusterFile_Tests(TestCase):
    """ Tests for UclustCreateClusterFile app controller """
    
    def setUp(self):
        
        self.tmp_sorted_fasta_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = "fasta")
        self.tmp_uc_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = "uc")
        
    def tearDown(self):
        if isfile(self.tmp_sorted_fasta_filepath):
            remove(self.tmp_sorted_fasta_filepath)
        if isfile(self.tmp_uc_filepath):
            remove(self.tmp_uc_filepath)
    
    def test_base_command(self):
        """ UclustCreateClusterFile should return the correct BaseCommand """
        c = UclustCreateClusterFile()
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "',getcwd(),'/"; ','uclust']))
        c.Parameters['--input'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "',getcwd(),'/"; ','uclust --input "seq.txt"']))
        c.Parameters['--uc'].on('sorted_output.fasta')
        # can't test against specific parameter order here, since 
        # params are written based on dict -- therefore just check
        # that changes to the parameter settings affect the base command
        # in the expected way
        c.Parameters['--id'].on(0.9)
        self.assertTrue('--id 0.9' in c.BaseCommand)
        c.Parameters['--rev'].on()
        self.assertTrue('--rev' in c.BaseCommand)
        c.Parameters['--rev'].off()
        self.assertFalse('--rev' in c.BaseCommand)

        

    def test_changing_working_dir(self):
        """ UclustCreateClusterFile BaseCommand should change dir
        
        Should change dir according to WorkingDir"""
        
        c = UclustCreateClusterFile(WorkingDir='/tmp/uclust_test')
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "','/tmp/uclust_test','/"; ','uclust']))
        c = UclustCreateClusterFile()
        c.WorkingDir = '/tmp/uclust_test2'
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "','/tmp/uclust_test2','/"; ','uclust']))
         
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/uclust_test')
        rmdir('/tmp/uclust_test2')
        
    def test_clusters_from_fasta_filepath(self):
        """ Should create clusters in uclust format from sorted fasta file 
        
        Since a fasta file has to be passed to the app controller for uclust,
        a temporary fasta file is created, and the sorted seqs supplied
        in this module are written to it.  This file is sent to the app 
        controller, and the resulting uclust file is compared to the expected
        results to ensure proper function of uclust as called by this app
        controller."""
        
        
        tmp_sorted_fasta = open(self.tmp_sorted_fasta_filepath,"w")
        for line in sorted_dna_seqs:
            tmp_sorted_fasta.write(line)

        tmp_sorted_fasta.close()

        test_app = UclustCreateClusterFile()
        test_app_res = test_app(data = \
         {'--input':self.tmp_sorted_fasta_filepath,\
         '--uc':self.tmp_uc_filepath, '--id':0.9})

        
        uc_file = open(test_app_res['ClusterFilepath'].name,"U")
        uc_file_res = []
        # Not appending comment lines of file, since the source data files
        # will change with each run, actual results are what we are 
        # interested in.
        for line in uc_file:
            if not(line.startswith("#")):
                parsed_line = "\t".join(line.split("\t")[:9])
                if not parsed_line.endswith('\n'):
                    parsed_line += '\n'
                uc_file_res.append(parsed_line)
            
        self.assertEqual(uc_file_res, uc_dna_clusters)
    
        test_app_res.cleanUp()
    
class UclustConvertToCdhit_Tests(TestCase):
    """ Tests for UclustConvertToCdhit app controller """
    
    def setUp(self):
        
        self.tmp_uc_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".uc")
        self.tmp_clstr_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = ".clstr")
        
    def tearDown(self):
        if isfile(self.tmp_uc_filepath):
            remove(self.tmp_uc_filepath)
        if isfile(self.tmp_clstr_filepath):
            remove(self.tmp_clstr_filepath)
    
    
    def test_base_command(self):
        """ UclustConvertToCdhit should return the correct BaseCommand """
        c = UclustConvertToCdhit()
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "',getcwd(),'/"; ','uclust']))
        c.Parameters['--uc2clstr'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "',getcwd(),'/"; ','uclust --uc2clstr "seq.txt"']))
        c.Parameters['--output'].on('sorted_output.clstr')
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "',getcwd(),'/"; ','uclust --uc2clstr "seq.txt" '+\
         '--output "sorted_output.clstr"']))

    def test_changing_working_dir(self):
        """ UclustConvertToCdhit BaseCommand should change WorkingDir"""
        
        c = UclustConvertToCdhit(WorkingDir='/tmp/uclust_test')
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "','/tmp/uclust_test','/"; ','uclust']))
        c = UclustConvertToCdhit()
        c.WorkingDir = '/tmp/uclust_test2'
        self.assertEqual(c.BaseCommand,\
         ''.join(['cd "','/tmp/uclust_test2','/"; ','uclust']))
         
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/uclust_test')
        rmdir('/tmp/uclust_test2')
        
    def test_convert_to_cdhit_from_uc_filepath(self):
        """ Should convert given uclust (.uc) file to cdhit (.clstr) format 
        
        Since a .uc file has to be passed to the app controller for uclust,
        a temporary .uc file is created, and the clusters supplied
        in this module are written to it.  This file is sent to the app 
        controller, and the resulting .clstr file is compared to the expected
        results to ensure proper function of uclust as called by this app
        controller."""
        

        tmp_uc = open(self.tmp_uc_filepath,"w")
        for line in uc_dna_clusters:
            # Need extra fields to be compatable with uclust 1.1
            tmp_uc.write(line.replace('\n','\t\n'))

        tmp_uc.close()

        test_app = UclustConvertToCdhit()
        
        
        
        test_app_res = test_app(data = \
           {'--uc2clstr':self.tmp_uc_filepath,'--output':self.tmp_clstr_filepath})

        
        clstr_file = open(test_app_res['CdhitFilepath'].name,"U")
        clstr_res = []
        for line in clstr_file:
            clstr_res.append(line.replace('\t',''))
            
        self.assertEqual(clstr_res, clstr_clusters)
   
        test_app_res.cleanUp()
    
class UclustSupporingModules(TestCase):
    """ Unit tests for supporting modules of uclust app controllers """

    def setUp(self):
        
        self.tmp_unsorted_fasta_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = "fasta")
        self.tmp_sorted_fasta_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = "fasta")
        self.tmp_uc_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = "uc")
        self.tmp_clstr_filepath = \
         get_tmp_filename(prefix = "uclust_test", suffix = "clstr")
        
        self.search_align_out1 = search_align_out1
        self.search_align_out1_expected = search_align_out1_expected
        self.search_align_query1_fp = \
         get_tmp_filename(prefix = "uclust_test", suffix = "clstr")
        open(self.search_align_query1_fp,'w').write(search_align_query1)
        self.search_align_template1_fp = \
         get_tmp_filename(prefix = "uclust_test", suffix = "clstr")
        open(self.search_align_template1_fp,'w').write(search_align_template1)
        
    def tearDown(self):
        if isfile(self.tmp_unsorted_fasta_filepath):
            remove(self.tmp_unsorted_fasta_filepath)
        if isfile(self.tmp_sorted_fasta_filepath):
            remove(self.tmp_sorted_fasta_filepath)
        if isfile(self.tmp_uc_filepath):
            remove(self.tmp_uc_filepath)
        remove(self.search_align_template1_fp)
        remove(self.search_align_query1_fp)
            

    def test_uclust_fasta_sort_from_filepath(self):
        """ Given an unsorted fasta filepath, will return sorted file """

        tmp_unsorted_fasta = open(self.tmp_unsorted_fasta_filepath,"w")
        for line in raw_dna_seqs:
            tmp_unsorted_fasta.write(line)

        tmp_unsorted_fasta.close()
        
        app_res = \
         uclust_fasta_sort_from_filepath(self.tmp_unsorted_fasta_filepath)
        
        sorted_fasta = open(app_res['SortedFasta'].name,"U")
        sorted_fasta_res = []
        for line in sorted_fasta:
            sorted_fasta_res.append(line)
            
        self.assertEqual(sorted_fasta_res, sorted_dna_seqs)
        app_res.cleanUp()
        
        
    def test_uclust_cluster_from_sorted_fasta_filepath(self):
        """ Given a sorted fasta filepath, will return uclust (.uc) file """
        
        tmp_sorted_fasta = open(self.tmp_sorted_fasta_filepath,"w")
        for line in sorted_dna_seqs:
            tmp_sorted_fasta.write(line)

        tmp_sorted_fasta.close()

        app_res = \
         uclust_cluster_from_sorted_fasta_filepath(self.tmp_sorted_fasta_filepath, \
         percent_ID = 0.90)

        
        uc_file = open(app_res['ClusterFilepath'].name,"U")
        uc_file_res = []
        # Not appending comment lines of file, since the source data files
        # will change with each run, actual results are what we are 
        # interested in.
        for line in uc_file:
            if not(line.startswith("#")):
                parsed_line = "\t".join(line.split("\t")[:9])
                if not parsed_line.endswith('\n'):
                    parsed_line += '\n'
                uc_file_res.append(parsed_line)
            
        self.assertEqual(uc_file_res, uc_dna_clusters)
        app_res.cleanUp()
        
    def test_uclust_convert_uc_to_cdhit_from_filepath(self):
        """ Given a uclust (.uc) file will return converted clstr file """
        
        tmp_uc = open(self.tmp_uc_filepath,"w")
        for line in uc_dna_clusters:
            # Need extra field to be compatable with uclust 1.1
            tmp_uc.write(line.replace('\n','\t\n'))

        tmp_uc.close()

        app_res = uclust_convert_uc_to_cdhit_from_filepath(self.tmp_uc_filepath)

        
        clstr_file = open(app_res['CdhitFilepath'].name,"U")
        clstr_res = []
        for line in clstr_file:
            clstr_res.append(line.replace('\t',''))
            
        self.assertEqual(clstr_res, clstr_clusters)
        app_res.cleanUp()
        
    def test_parse_uclust_clstr_file(self):
        """ Ensures that list of lists of OTUs will be returned """
        
        clusters_res = parse_uclust_clstr_file(clstr_clusters)
        
        self.assertEqual(clusters_res, expected_cluster_list)
    
    def test_get_output_filepaths(self):
        """ Properly generates output filepath names """
        
        fasta_res, uc_res, cd_hit_res, output_dir_res = \
         get_output_filepaths("/tmp/","test_seqs.fasta")
        
        self.assertEqual(fasta_res, "/tmp/test_seqs_sorted.fasta")
        self.assertEqual(uc_res, "/tmp/test_seqs_sorted.uc")
        self.assertEqual(cd_hit_res, "/tmp/test_seqs_cdhit.clstr")
        self.assertEqual(output_dir_res, "/tmp")
        
        fasta_res, uc_res, cd_hit_res, output_dir_res = \
         get_output_filepaths(".","test_seqs.fasta")
        self.assertEqual(fasta_res, ".//test_seqs_sorted.fasta")
        self.assertEqual(uc_res, ".//test_seqs_sorted.uc")
        self.assertEqual(cd_hit_res, ".//test_seqs_cdhit.clstr")
        self.assertEqual(output_dir_res, "./")
        
    def test_get_clusters_from_fasta_filepath(self):
        """ Tests for return of lists of OTUs from given fasta filepath """
        
        tmp_unsorted_fasta = open(self.tmp_unsorted_fasta_filepath,"w")
        for line in raw_dna_seqs:
            tmp_unsorted_fasta.write(line)

        tmp_unsorted_fasta.close()
        
        clusters_res = \
         get_clusters_from_fasta_filepath(self.tmp_unsorted_fasta_filepath, \
          percent_ID = 0.90)

        self.assertEqual(clusters_res, expected_cluster_list)
        
    def test_process_uclust_blast_result(self):
        """parsing of pairwise alignments functions as expected """
        actual = list(process_uclust_blast_result(self.search_align_out1))
        expected = self.search_align_out1_expected
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

raw_dna_seqs = ['>uclust_test_seqs_0\n',
'ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT\n',
'>uclust_test_seqs_1\n',
'GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT\n',
'>uclust_test_seqs_2\n',
'CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT\n',
'>uclust_test_seqs_3\n',
'CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT\n',
'>uclust_test_seqs_4\n',
'ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT\n',
'>uclust_test_seqs_5\n',
'CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA\n',
'>uclust_test_seqs_6\n',
'CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT\n',
'>uclust_test_seqs_7\n',
'AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT\n',
'>uclust_test_seqs_8\n',
'CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT\n',
'>uclust_test_seqs_9\n',
'GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA']

sorted_dna_seqs=['>uclust_test_seqs_7\n',
'AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT\n',
'>uclust_test_seqs_4\n',
'ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT\n',
'>uclust_test_seqs_2\n',
'CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT\n',
'>uclust_test_seqs_3\n',
'CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT\n',
'>uclust_test_seqs_1\n',
'GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT\n',
'>uclust_test_seqs_5\n',
'CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA\n',
'>uclust_test_seqs_6\n',
'CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT\n',
'>uclust_test_seqs_0\n',
'ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT\n',
'>uclust_test_seqs_8\n',
'CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT\n',
'>uclust_test_seqs_9\n',
'GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA\n']

# Clusters are created at a 0.90% identity
uc_dna_clusters=[
'S	0	80	*	*	*	*	*	uclust_test_seqs_7\n',
'S	1	79	*	*	*	*	*	uclust_test_seqs_4\n',
'S	2	78	*	*	*	*	*	uclust_test_seqs_2\n',
'S	3	77	*	*	*	*	*	uclust_test_seqs_3\n',
'S	4	76	*	*	*	*	*	uclust_test_seqs_1\n',
'S	5	75	*	*	*	*	*	uclust_test_seqs_5\n',
'S	6	74	*	*	*	*	*	uclust_test_seqs_6\n',
'S	7	73	*	*	*	*	*	uclust_test_seqs_0\n',
'H	6	72	91.7	+	0	0	2I72M	uclust_test_seqs_8\n',
'S	8	71	*	*	*	*	*	uclust_test_seqs_9\n',
'C	0	1	*	*	*	*	*	uclust_test_seqs_7\n',
'C	1	1	*	*	*	*	*	uclust_test_seqs_4\n',
'C	2	1	*	*	*	*	*	uclust_test_seqs_2\n',
'C	3	1	*	*	*	*	*	uclust_test_seqs_3\n',
'C	4	1	*	*	*	*	*	uclust_test_seqs_1\n',
'C	5	1	*	*	*	*	*	uclust_test_seqs_5\n',
'C	6	2	91.7	*	*	*	*	uclust_test_seqs_6\n',
'C	7	1	*	*	*	*	*	uclust_test_seqs_0\n',
'C	8	1	*	*	*	*	*	uclust_test_seqs_9\n']

"""
# Old version, incompatible with uclust 1.1 format
# Clusters are created at a 0.90% identity
uc_dna_clusters=[
'S	0	80	*	*	*	*	*	uclust_test_seqs_7	*\n',
'S	1	79	*	*	*	*	*	uclust_test_seqs_4	*\n',
'S	2	78	*	*	*	*	*	uclust_test_seqs_2	*\n',
'S	3	77	*	*	*	*	*	uclust_test_seqs_3	*\n',
'S	4	76	*	*	*	*	*	uclust_test_seqs_1	*\n',
'S	5	75	*	*	*	*	*	uclust_test_seqs_5	*\n',
'S	6	74	*	*	*	*	*	uclust_test_seqs_6	*\n',
'S	7	73	*	*	*	*	*	uclust_test_seqs_0	*\n',
'H	6	72	91.7	+	0	0	2I72M	uclust_test_seqs_8	*\n',
'S	8	71	*	*	*	*	*	uclust_test_seqs_9	*\n',
'C	0	1	*	*	*	*	*	uclust_test_seqs_7	*\n',
'C	1	1	*	*	*	*	*	uclust_test_seqs_4	*\n',
'C	2	1	*	*	*	*	*	uclust_test_seqs_2	*\n',
'C	3	1	*	*	*	*	*	uclust_test_seqs_3	*\n',
'C	4	1	*	*	*	*	*	uclust_test_seqs_1	*\n',
'C	5	1	*	*	*	*	*	uclust_test_seqs_5	*\n',
'C	6	2	91.7	*	*	*	*	uclust_test_seqs_6	*\n',
'C	7	1	*	*	*	*	*	uclust_test_seqs_0	*\n',
'C	8	1	*	*	*	*	*	uclust_test_seqs_9	*\n'] """

clstr_clusters=['>Cluster 0\n',
'0       80nt, >uclust_test_seqs_7... *\n',
'>Cluster 1\n',
'0       79nt, >uclust_test_seqs_4... *\n',
'>Cluster 2\n',
'0       78nt, >uclust_test_seqs_2... *\n',
'>Cluster 3\n',
'0       77nt, >uclust_test_seqs_3... *\n',
'>Cluster 4\n',
'0       76nt, >uclust_test_seqs_1... *\n',
'>Cluster 5\n',
'0       75nt, >uclust_test_seqs_5... *\n',
'>Cluster 6\n',
'0       74nt, >uclust_test_seqs_6... *\n',
'1       72nt, >uclust_test_seqs_8... at +/92%\n',
'>Cluster 7\n',
'0       73nt, >uclust_test_seqs_0... *\n',
'>Cluster 8\n',
'0       71nt, >uclust_test_seqs_9... *\n']

expected_cluster_list=[['uclust_test_seqs_7'], ['uclust_test_seqs_4'], ['uclust_test_seqs_2'], ['uclust_test_seqs_3'], ['uclust_test_seqs_1'], ['uclust_test_seqs_5'], ['uclust_test_seqs_6', 'uclust_test_seqs_8'], ['uclust_test_seqs_0'], ['uclust_test_seqs_9']]

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

search_align_out1 = """# uclust --input sm_query.fasta --lib sm_template.fasta --id 0.75 --libonly --rev --maxaccepts 0 --blastout s_bl_out.txt --blast_termgaps
# version=1.1.572

Query  >1_like
Target >1

  1 + -------------------------------TACGGCTACCTTGTTACGACTTCATCCCAATCA 33
                                     |||||||||||||||||||||||| ||||||||
  1 + AGAAAGGAGGTGATCCAGCCGCACCTTCCGATACGGCTACCTTGTTACGACTTCACCCCAATCA 64

 34 + TTTGTTCCACCTTCGACGGCTA------------------------------------------ 55
      ||||||||||||||||||||||                                          
 65 + TTTGTTCCACCTTCGACGGCTAGCTCCAAATGGTTACTCCACCGGCTTCGGGTGTTACAAACTC 128

Identities 54/55 (98.2%), gaps 73/128 (57.0%), Id 54/55 (98.2%)

Query  >2_like
Target >2

 1 + -------------------ATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGG 45
                        |||||||||||||||||||||||||||||||||||||||||||||
 1 + AGCCCAAATCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGG 64

46 + GATGGCAACTAAG--------------- 58
     |||||||||||||               
65 + GATGGCAACTAAGCTTAAGGGTTGCGCT 92

Identities 58/58 (100.0%), gaps 34/92 (37.0%), Id 58/58 (100.0%)

Query  >2_like_rc
Target >2

 1 + ---------------CTTAGTTGCCATCCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCA 49
                    |||||||||||||||||||||||||||||||||||||||||||||||||
92 - AGCGCAACCCTTAAGCTTAGTTGCCATCCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCA 29

50 + AATCATCAT------------------- 58
     |||||||||                   
28 - AATCATCATGCCCCTTATGATTTGGGCT 1

Identities 58/58 (100.0%), gaps 34/92 (37.0%), Id 58/58 (100.0%)
""".split('\n')

search_align_out1_expected = [
         ('1_like','1','-------------------------------TACGGCTACCTTGTTACGACTTCATCCCAATCATTTGTTCCACCTTCGACGGCTA------------------------------------------','AGAAAGGAGGTGATCCAGCCGCACCTTCCGATACGGCTACCTTGTTACGACTTCACCCCAATCATTTGTTCCACCTTCGACGGCTAGCTCCAAATGGTTACTCCACCGGCTTCGGGTGTTACAAACTC',98.2),
         
         ('2_like','2','-------------------ATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAG---------------','AGCCCAAATCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAGCTTAAGGGTTGCGCT',100.0),\
         
         ('2_like_rc RC','2','-------------------ATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAG---------------','AGCCCAAATCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGGATGGCAACTAAGCTTAAGGGTTGCGCT',100.0)]
         
if __name__ == '__main__':
    main()
