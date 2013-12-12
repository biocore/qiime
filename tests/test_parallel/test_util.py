#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import exists
from cogent.util.unit_test import TestCase, main
from cogent import LoadSeqs
from cogent.util.misc import remove_files
from qiime.util import (get_tmp_filename,
                        get_qiime_temp_dir) 
from qiime.parallel.util import (ParallelWrapper,
                                 BufferedWriter)

class ParallelWrapperTests(TestCase):

    def setUp(self):
        """ """
        # instantiating the abstract class to test some of the more 
        # stand-alone methods
        self.pw = ParallelWrapper()

    def test_merge_to_n_commands_even(self):
        """ _merge_to_n_commands functions as expected (even number of cmds)"""
        commands = ['pick_otus.py -h ; mv somthing.txt something_else.txt',
                    'pick_otus.py -g',
                    'pick_otus.py -f',
                    'pick_otus.py -w']
                    
        expected = ['/bin/bash ; pick_otus.py -h ; mv somthing.txt something_else.txt ; pick_otus.py -g ; pick_otus.py -f ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,1)
        self.assertEqual(actual,expected)
        
        expected = [
         '/bin/bash ; pick_otus.py -h ; mv somthing.txt something_else.txt ; pick_otus.py -g ; exit',
         '/bin/bash ; pick_otus.py -f ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,2)
        self.assertEqual(actual,expected)
        
        # rounds to 2 jobs to start
        expected = [
         '/bin/bash ; pick_otus.py -h ; mv somthing.txt something_else.txt ; pick_otus.py -g ; exit',
         '/bin/bash ; pick_otus.py -f ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,3)
        self.assertEqual(actual,expected)
        
        expected = ['/bin/bash ; pick_otus.py -h ; mv somthing.txt something_else.txt ; exit',
                    '/bin/bash ; pick_otus.py -g ; exit',
                    '/bin/bash ; pick_otus.py -f ; exit',
                    '/bin/bash ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,4)
        self.assertEqual(actual,expected)
        
        self.assertRaises(ValueError,self.pw._merge_to_n_commands,commands,0)
        self.assertRaises(ValueError,self.pw._merge_to_n_commands,commands,-42)
        
        # jobs to start is much higer than actual jobs
        expected = ['/bin/bash ; pick_otus.py -h ; mv somthing.txt something_else.txt ; exit',
                    '/bin/bash ; pick_otus.py -g ; exit',
                    '/bin/bash ; pick_otus.py -f ; exit',
                    '/bin/bash ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,100)
        self.assertEqual(actual,expected)
        
        self.assertRaises(ValueError,self.pw._merge_to_n_commands,commands,0)
        self.assertRaises(ValueError,self.pw._merge_to_n_commands,commands,-42)
        
        
    def test_merge_to_n_commands_odd(self):
        """ _merge_to_n_commands functions as expected (odd number of cmds)"""
        commands = ['pick_otus.py -h',
                    'pick_otus.py -g',
                    'pick_otus.py -w']
                    
        expected = ['/bin/bash ; pick_otus.py -h ; pick_otus.py -g ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,1)
        self.assertEqual(actual,expected)
                    
        # rounds to 1 job to start
        expected = ['/bin/bash ; pick_otus.py -h ; pick_otus.py -g ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,2)
        self.assertEqual(actual,expected)
                    
        expected = ['/bin/bash ; pick_otus.py -h ; exit',
                    '/bin/bash ; pick_otus.py -g ; exit',
                    '/bin/bash ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,3)
        self.assertEqual(actual,expected)
        
        expected = ['/bin/bash ; pick_otus.py -h ; exit',
                    '/bin/bash ; pick_otus.py -g ; exit',
                    '/bin/bash ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,4)
        self.assertEqual(actual,expected)
        
        expected = ['/bin/bash ; pick_otus.py -h ; exit',
                    '/bin/bash ; pick_otus.py -g ; exit',
                    '/bin/bash ; pick_otus.py -w ; exit']
        actual = self.pw._merge_to_n_commands(commands,100)
        self.assertEqual(actual,expected)
        
        self.assertRaises(ValueError,self.pw._merge_to_n_commands,commands,0)
        self.assertRaises(ValueError,self.pw._merge_to_n_commands,commands,-42)    
        
    def test_merge_to_n_commands_alt_params(self):
        """ _merge_to_n_commands functions with alt params"""
        commands = ['pick_otus.py -h',
                    'pick_otus.py -g',
                    'pick_otus.py -w']
                    
        expected = ['pick_otus.py -h ; pick_otus.py -g ; pick_otus.py -w']
        actual = self.pw._merge_to_n_commands(commands,2,command_prefix='',command_suffix='')
        self.assertEqual(actual,expected)
                    
        expected = ['pick_otus.py -h ! pick_otus.py -g ! pick_otus.py -w']
        actual = self.pw._merge_to_n_commands(commands,2,command_prefix='',
         command_suffix='',delimiter=' ! ')
        self.assertEqual(actual,expected)
        
        commands = map(str,range(10))
        actual = self.pw._merge_to_n_commands(commands,5,command_prefix='',
         command_suffix='',delimiter=',')
        expected = ['0,1','2,3','4,5','6,7','8,9']
        self.assertEqual(actual,expected)

    def test_merge_to_n_commands_w_prefix(self):
        """ _merge_to_n_commands functions as expected (w prefix/suffix)"""
        commands = ['/bin/bash ; pick_otus.py -h ; exit',
                    '/bin/bash;pick_otus.py -g;exit',
                    '/bin/bash ; pick_otus.py -h ; exit ; /bin/bash ; pick_otus.py -w ; exit']
        expected = ['/bin/bash ; pick_otus.py -h ; pick_otus.py -g ; pick_otus.py -h ; pick_otus.py -w ; exit']
        
        actual = self.pw._merge_to_n_commands(commands,2,command_prefix='/bin/bash ;',command_suffix='; exit')
        self.assertEqual(actual,expected)
        actual = self.pw._merge_to_n_commands(commands,2)
        self.assertEqual(actual,expected)

    def test_get_random_job_prefix(self):
        """ get_random_job_prefix functions as expected """
        
        s1 = self.pw._get_random_job_prefix()
        s2 = self.pw._get_random_job_prefix()
        self.assertNotEqual(s1,s2)
        self.assertEqual(len(s1),10)
        self.assertEqual(len(s2),10)
        
        # different max len
        s1 = self.pw._get_random_job_prefix(max_job_prefix_len=22)
        self.assertEqual(len(s1),22)
        
        # fixed_prefix added
        s1 = self.pw._get_random_job_prefix(fixed_prefix='TEST')
        s2 = self.pw._get_random_job_prefix(fixed_prefix='TEST')
        self.assertNotEqual(s1,s2)
        self.assertEqual(len(s1),10)
        self.assertTrue(s1.startswith('TEST'))
        self.assertTrue(s2.startswith('TEST'))
        # leading/trailing underscores added
        self.assertTrue(s1.startswith('TEST_'))
        self.assertTrue(s1.endswith('_'))
        
        # no leading/trailing underscores
        s1 = self.pw._get_random_job_prefix(leading_trailing_underscores=False)
        self.assertFalse(s1.startswith('_'))
        self.assertFalse(s1.endswith('_'))
        
        # combo of all parameters
        s1 = self.pw._get_random_job_prefix(leading_trailing_underscores=False,\
         fixed_prefix='HELLO',max_job_prefix_len=12)
        self.assertEqual(len(s1),12)
        self.assertTrue(s1.startswith('HELLO'))
        self.assertFalse(s1.endswith('_'))
    
    def test_compute_seqs_per_file(self):
        """compute_seqs_per_file functions as expected
        """
        temp_fasta_fp = get_tmp_filename(\
         prefix='QiimeScriptUtilTests',suffix='.fasta')
        temp_fasta = ['>seq','AAACCCCAAATTGG'] * 25
        open(temp_fasta_fp,'w').write('\n'.join(temp_fasta))
        
        actual_25 = self.pw._compute_seqs_per_file(temp_fasta_fp,25)
        actual_2 = self.pw._compute_seqs_per_file(temp_fasta_fp,2)
        actual_10 = self.pw._compute_seqs_per_file(temp_fasta_fp,10)
        actual_5 = self.pw._compute_seqs_per_file(temp_fasta_fp,5)
        actual_40 = self.pw._compute_seqs_per_file(temp_fasta_fp,40)
        
        remove_files([temp_fasta_fp])
        
        self.assertEqual(actual_25,1)
        self.assertEqual(actual_2,13)
        self.assertEqual(actual_10,3)
        self.assertEqual(actual_5,5)
        self.assertEqual(actual_40,1)

class BufferedWriterTests(TestCase):
    
    def setUp(self):
        """ """
        self.files_to_remove = []
        tmp_dir = get_qiime_temp_dir()
        self.test_fp = get_tmp_filename(tmp_dir=tmp_dir,
                                        prefix='bufWriterTest',
                                        suffix='.txt')
        self.files_to_remove.append(self.test_fp)

    def tearDown(self):
        """ """
        remove_files(self.files_to_remove)

    def test_init(self):
        """BufferedWriter constructor works"""
        
        b = BufferedWriter(self.test_fp)
        self.assertTrue(exists(self.test_fp))

    def test_write(self):
        """BufferedWriter writes nothing until max buffer reached."""
        
        b = BufferedWriter(self.test_fp, buf_size=2)
        b.write("1")
        content = open(self.test_fp, "r").readlines()
        self.assertEquals(content, [])

        #still nothing
        b.write("2")
        content = open(self.test_fp, "r").readlines()
        self.assertEquals(content, [])

        #finally, buffer is flushed
        b.write("3")
        content = open(self.test_fp, "r").readlines()
        self.assertEquals(content, ["123"])

    def test_close(self):
        """close() flushes the buffer"""
        
        b = BufferedWriter(self.test_fp, buf_size=2)
        b.write("1")
        b.close()
        content = open(self.test_fp, "r").readlines()
        self.assertEquals(content, ["1"])

if __name__ == "__main__":
    main()
