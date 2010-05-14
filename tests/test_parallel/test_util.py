#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# test_util.py

""" Description
File created on 25 Aug 2009.

"""
from __future__ import division
from os import remove
from cogent import LoadSeqs
from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
from qiime.parallel.util import split_fasta, get_random_job_prefix,\
 write_jobs_file, compute_seqs_per_file, build_filepaths_from_filepaths,\
 submit_jobs

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

class UtilTests(TestCase):
    """Tests of parallel code utility functions """
    
    def setUp(self):
        pass
        
    def test_get_random_job_prefix(self):
        """ get_random_job_prefix functions as expected """
        
        s1 = get_random_job_prefix()
        s2 = get_random_job_prefix()
        self.assertNotEqual(s1,s2)
        self.assertEqual(len(s1),10)
        self.assertEqual(len(s2),10)
        
        # different max len
        s1 = get_random_job_prefix(max_job_prefix_len=22)
        self.assertEqual(len(s1),22)
        
        # fixed_prefix added
        s1 = get_random_job_prefix(fixed_prefix='TEST')
        s2 = get_random_job_prefix(fixed_prefix='TEST')
        self.assertNotEqual(s1,s2)
        self.assertEqual(len(s1),10)
        self.assertTrue(s1.startswith('TEST'))
        self.assertTrue(s2.startswith('TEST'))
        # leading/trailing underscores added
        self.assertTrue(s1.startswith('TEST_'))
        self.assertTrue(s1.endswith('_'))
        
        # no leading/trailing underscores
        s1 = get_random_job_prefix(leading_trailing_underscores=False)
        self.assertFalse(s1.startswith('_'))
        self.assertFalse(s1.endswith('_'))
        
        # combo of all parameters
        s1 = get_random_job_prefix(leading_trailing_underscores=False,\
         fixed_prefix='HELLO',max_job_prefix_len=12)
        self.assertEqual(len(s1),12)
        self.assertTrue(s1.startswith('HELLO'))
        self.assertFalse(s1.endswith('_'))
        
    def test_submit_jobs_fail(self):
        """submit jobs fails by raising an error
        """
        self.assertRaises(RuntimeError,submit_jobs,
         'some_fake_exe','some_fake_fp.txt','JOB')
        
    def test_compute_seqs_per_file(self):
        """compute_seqs_per_file functions as expected
        """
        temp_fasta_fp = get_tmp_filename(\
         prefix='QiimeScriptUtilTests',suffix='.fasta')
        temp_fasta = ['>seq','AAACCCCAAATTGG'] * 25
        open(temp_fasta_fp,'w').write('\n'.join(temp_fasta))
        
        actual_25 = compute_seqs_per_file(temp_fasta_fp,25)
        actual_2 = compute_seqs_per_file(temp_fasta_fp,2)
        actual_10 = compute_seqs_per_file(temp_fasta_fp,10)
        actual_5 = compute_seqs_per_file(temp_fasta_fp,5)
        actual_40 = compute_seqs_per_file(temp_fasta_fp,40)
        
        remove(temp_fasta_fp)
        
        self.assertEqual(actual_25,1)
        self.assertEqual(actual_2,13)
        self.assertEqual(actual_10,3)
        self.assertEqual(actual_5,5)
        self.assertEqual(actual_40,1)
        
        
    def test_build_filepaths_from_filepaths(self):
        """ build_filepaths_from_filepaths functions as expected
        """
        # no additional params just strips paths
        in_fps = ['in1.txt','somewhere/in2.txt','in3.txt']
        actual = build_filepaths_from_filepaths(in_fps)
        expected = ['in1.txt','in2.txt','in3.txt']
        self.assertEqual(actual,expected)
        
        # replace works as expected
        in_fps = ['in1.txt','somewhere/in2.txt','in3.txt']
        actual = build_filepaths_from_filepaths(in_fps,replacement=('.txt','.fasta'))
        expected = ['in1.fasta','in2.fasta','in3.fasta']
        self.assertEqual(actual,expected)
        
        # adding directory works as expected (no trailing / in directory)
        in_fps = ['in1.txt','somewhere/in2.txt','in3.txt']
        actual = build_filepaths_from_filepaths(in_fps,directory='/home/bob')
        expected = ['/home/bob/in1.txt','/home/bob/in2.txt','/home/bob/in3.txt']
        self.assertEqual(actual,expected)
        # adding directory works as expected (trailing / in directory)
        in_fps = ['in1.txt','somewhere/in2.txt','in3.txt']
        actual = build_filepaths_from_filepaths(in_fps,directory='/home/bob/')
        expected = ['/home/bob/in1.txt','/home/bob/in2.txt','/home/bob/in3.txt']
        self.assertEqual(actual,expected)
        
        # prefix works as expected
        in_fps = ['in1.txt','somewhere/in2.txt','in3.txt']
        actual = build_filepaths_from_filepaths(in_fps,prefix='blah_')
        expected = ['blah_in1.txt','blah_in2.txt','blah_in3.txt']
        self.assertEqual(actual,expected)
        
        # suffix works as expected
        in_fps = ['in1.txt','somewhere/in2.txt','in3.txt']
        actual = build_filepaths_from_filepaths(in_fps,suffix='.out')
        expected = ['in1.txt.out','in2.txt.out','in3.txt.out']
        self.assertEqual(actual,expected)
        
        # combination works as expected (including order of operations)
        in_fps = ['txt_dir/in1.txt','somewhere/in2.txt','txt_dir/in3.txt']
        actual = build_filepaths_from_filepaths(in_fps,replacement=('txt','fasta'),\
         prefix='blah_',suffix='.out',directory='out/')
        expected = ['out/blah_in1.fasta.out','out/blah_in2.fasta.out',\
                    'out/blah_in3.fasta.out']
        self.assertEqual(actual,expected)
    
        actual = build_filepaths_from_filepaths(in_fps,prefix='',\
            directory='',suffix='',replacement=('',''))
    
    def test_split_fasta_equal_num_seqs_per_file(self):
        """split_fasta funcs as expected when equal num seqs go to each file
        """
        filename_prefix = get_random_job_prefix(fixed_prefix='/tmp/')
        infile = ['>seq1','AACCTTAA','>seq2','TTAACC','AATTAA',\
         '>seq3','CCTT--AA']
         
        actual = split_fasta(infile, 1, filename_prefix)
        actual_seqs = []
        for fp in actual:
            actual_seqs += list(open(fp))
        remove_files(actual)
        
        expected = ['%s.%d.fasta' % (filename_prefix,i) for i in range(3)]
        
        self.assertEqual(actual,expected)
        self.assertEqual(\
         LoadSeqs(data=infile,aligned=False),\
         LoadSeqs(data=actual_seqs,aligned=False))
        
        
    def test_split_fasta_diff_num_seqs_per_file(self):
        """split_fasta funcs as expected when diff num seqs go to each file
        """
        filename_prefix = get_random_job_prefix(fixed_prefix='/tmp/')
        infile = ['>seq1','AACCTTAA','>seq2','TTAACC','AATTAA',\
         '>seq3','CCTT--AA']
         
        actual = split_fasta(infile, 2, filename_prefix)
        
        actual_seqs = []
        for fp in actual:
            actual_seqs += list(open(fp))
        remove_files(actual)
        
        expected = ['%s.%d.fasta' % (filename_prefix,i) for i in range(2)]
        # list of file paths is as expected
        self.assertEqual(actual,expected)
        # building seq collections from infile and the split files result in
        # equivalent seq collections
        self.assertEqual(\
         LoadSeqs(data=infile,aligned=False),\
         LoadSeqs(data=actual_seqs,aligned=False))
         
    def test_split_fasta_diff_num_seqs_per_file_alt(self):
        """split_fasta funcs always catches all seqs
        """
        # start with 59 seqs (b/c it's prime, so should make more 
        # confusing splits)
        in_seqs = LoadSeqs(data=[('seq%s' % k,'AACCTTAA') for k in range(59)])
        infile = in_seqs.toFasta().split('\n')
        
        # test seqs_per_file from 1 to 1000
        for i in range(1,1000):
            filename_prefix = get_random_job_prefix(fixed_prefix='/tmp/')
         
            actual = split_fasta(infile, i, filename_prefix)
        
            actual_seqs = []
            for fp in actual:
                actual_seqs += list(open(fp))
            # remove the files now, so if the test fails they still get 
            # cleaned up
            remove_files(actual)
            
            # building seq collections from infile and the split files result in
            # equivalent seq collections
            self.assertEqual(\
             LoadSeqs(data=infile,aligned=False),\
             LoadSeqs(data=actual_seqs,aligned=False))
         
    
        
if __name__ == "__main__":
    main()
