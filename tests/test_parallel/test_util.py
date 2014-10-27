#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import close
from os.path import exists
from shutil import rmtree
from tempfile import mkstemp, mkdtemp
from unittest import TestCase, main

from skbio.util.misc import remove_files

from qiime.util import get_qiime_temp_dir
from qiime.parallel.util import (BufferedWriter, concatenate_files,
                                 merge_files_from_dirs, input_fasta_splitter)


class UtilTests(TestCase):
    def setUp(self):
        self.files_to_remove = []
        self.dirs_to_remove = []

        self.tmp_dir = get_qiime_temp_dir()

        # Create some files to concatenate them
        fd, self.fp1 = mkstemp(dir=self.tmp_dir)
        close(fd)
        with open(self.fp1, 'w') as f:
            f.write(FILE1)
        self.files_to_remove.append(self.fp1)

        fd, self.fp2 = mkstemp(dir=self.tmp_dir)
        close(fd)
        with open(self.fp2, 'w') as f:
            f.write(FILE2)
        self.files_to_remove.append(self.fp2)

        fd, self.fp3 = mkstemp(dir=self.tmp_dir)
        close(fd)
        with open(self.fp3, 'w') as f:
            f.write(FILE3)
        self.files_to_remove.append(self.fp3)

        # Create some directories with some files so we can test
        # merge files from dirs
        self.dir1 = mkdtemp(dir=self.tmp_dir)
        self.dirs_to_remove.append(self.dir1)
        fd, fp = mkstemp(dir=self.dir1, suffix="_test.txt")
        close(fd)
        with open(fp, 'w') as f:
            f.write(FILE1)
        fd, fp = mkstemp(dir=self.dir1, suffix="_test.log")
        close(fd)
        with open(fp, 'w') as f:
            f.write(LOGFILE1)

        self.dir2 = mkdtemp(dir=self.tmp_dir)
        self.dirs_to_remove.append(self.dir2)
        fd, fp = mkstemp(dir=self.dir2, suffix="_test.txt")
        close(fd)
        with open(fp, 'w') as f:
            f.write(FILE2)
        fd, fp = mkstemp(dir=self.dir2, suffix="_test.log")
        close(fd)
        with open(fp, 'w') as f:
            f.write(LOGFILE2)

        self.dir3 = mkdtemp(dir=self.tmp_dir)
        self.dirs_to_remove.append(self.dir3)
        fd, fp = mkstemp(dir=self.dir3, suffix="_test.txt")
        close(fd)
        with open(fp, 'w') as f:
            f.write(FILE3)
        fd, fp = mkstemp(dir=self.dir3, suffix="_test.log")
        close(fd)
        with open(fp, 'w') as f:
            f.write(LOGFILE3)

        # Create a fasta file with some sequences on it so we can test
        # the input fasta splitter
        fd, self.fasta_fp = mkstemp(dir=self.tmp_dir, prefix="first_",
                                    suffix=".fasta")
        self.files_to_remove.append(self.fasta_fp)
        with open(self.fasta_fp, 'w') as f:
            f.write(FASTA_FILE)

        # This file is also for testing the input_fasta_splitter, in the
        # previous version of the input_fasta_splitter, this test will return
        # 3 files instead of 4, making the current parallel framework to crash
        fd, self.fasta_fp_2 = mkstemp(dir=self.tmp_dir, prefix="second_",
                                      suffix=".fasta")
        self.files_to_remove.append(self.fasta_fp_2)
        with open(self.fasta_fp_2, 'w') as f:
            f.write(FASTA_FILE2)

    def tearDown(self):
        remove_files(self.files_to_remove)
        # Remove directories last, so we don't get errors trying to remove
        # files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_concatenate_files(self):
        fd, out_file = mkstemp(dir=self.tmp_dir)
        close(fd)
        self.files_to_remove.append(out_file)

        concatenate_files(out_file, [self.fp1, self.fp2, self.fp3])

        with open(out_file, 'U') as f:
            obs = f.read()

        self.assertEqual(obs, CONCAT_RESULT)

    def test_merge_files_from_dirs(self):
        fd, out_file = mkstemp(dir=self.tmp_dir)
        close(fd)
        self.files_to_remove.append(out_file)

        merge_files_from_dirs(out_file, [self.dir1, self.dir2, self.dir3],
                              "*_test.txt", concatenate_files)

        with open(out_file, 'U') as f:
            obs = f.read()

        self.assertEqual(obs, CONCAT_RESULT)

    def test_input_fasta_splitter_even(self):
        out_dir = mkdtemp(dir=self.tmp_dir)
        self.dirs_to_remove.append(out_dir)

        paths = input_fasta_splitter(self.fasta_fp, out_dir, 2)

        # Check that the number of paths returned is correct
        self.assertEqual(len(paths), 2)

        # Check that the paths exists
        for p in paths:
            self.assertTrue(exists(p))

        # Check the contents of each file
        with open(paths[0], 'U') as f:
            obs = f.read()
        exp = """>11472286
GATGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGAACACTTTGTGTTTTGAGTTAATAGTTCGATAGTA
>11468680
TAAACTGAAGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACG
>11469739
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGAGAAGCTAACTTCTGA
>11460523
AGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGCGAAATCGGGCACTCA
>11480235
TGGTTTGATCCTGGCTCAGGATTAACGCTGGCGGCGCGCCTTATACATGCAAGTCGAACGAGCCTTGTGCTTCGCAC
"""
        self.assertEqual(obs, exp)

        with open(paths[1], 'U') as f:
            obs = f.read()
        exp = """>11472384
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACGGGGGCAAC
>11458037
GACGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGGTTTCGAAGATCGGACTTCGAATTTCGAATTTCGAT
>11469752
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGCAGCGAGTTCCTCAC
>11460543
TGGTTTGATCCTGGCTCAGGACAAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGAGAAGCCAGCTTTTGAT
>11480408
AATTTAGCGGCCGCGAATTCGCCCTTGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCA
"""
        self.assertEqual(obs, exp)

    def test_input_fasta_splitter_uneven(self):
        out_dir = mkdtemp(dir=self.tmp_dir)
        self.dirs_to_remove.append(out_dir)

        paths = input_fasta_splitter(self.fasta_fp, out_dir, 3)

        # Check that the number of paths returned is correct
        self.assertEqual(len(paths), 3)

        # Check that the paths exists
        for p in paths:
            self.assertTrue(exists(p))

        # Check the contents of each file
        with open(paths[0], 'U') as f:
            obs = f.read()
        exp = """>11472286
GATGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGAACACTTTGTGTTTTGAGTTAATAGTTCGATAGTA
>11458037
GACGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGGTTTCGAAGATCGGACTTCGAATTTCGAATTTCGAT
>11460523
AGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGCGAAATCGGGCACTCA
>11480408
AATTTAGCGGCCGCGAATTCGCCCTTGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCA
"""
        self.assertEqual(obs, exp)

        with open(paths[1], 'U') as f:
            obs = f.read()
        exp = """>11472384
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACGGGGGCAAC
>11469739
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGAGAAGCTAACTTCTGA
>11460543
TGGTTTGATCCTGGCTCAGGACAAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGAGAAGCCAGCTTTTGAT
"""
        self.assertEqual(obs, exp)

        with open(paths[2], 'U') as f:
            obs = f.read()
        exp = """>11468680
TAAACTGAAGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACG
>11469752
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGCAGCGAGTTCCTCAC
>11480235
TGGTTTGATCCTGGCTCAGGATTAACGCTGGCGGCGCGCCTTATACATGCAAGTCGAACGAGCCTTGTGCTTCGCAC
"""
        self.assertEqual(obs, exp)

    def test_input_fasta_splitter_special(self):
        # This test will fail with the previous implementation of the input
        # fasta splitter, which will break the current parallel implementation
        # Since it is a special case, we add a specific test for it, so further
        # changes in the code does not break the parallel framework
        out_dir = mkdtemp(dir=self.tmp_dir)
        self.dirs_to_remove.append(out_dir)

        paths = input_fasta_splitter(self.fasta_fp_2, out_dir, 4)

        # Check that the number of paths returned is correct
        self.assertEqual(len(paths), 4)

        # Check that the paths exists
        for p in paths:
            self.assertTrue(exists(p))

        # check the contents of each file

        with open(paths[0], 'U') as f:
            obs = f.read()
        exp = """>11472286
GATGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGAACACTTTGTGTTTTGAGTTAATAGTTCGATAGTA
>11469739
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGAGAAGCTAACTTCTGA
"""
        self.assertEqual(obs, exp)

        with open(paths[1], 'U') as f:
            obs = f.read()
        exp = """>11472384
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACGGGGGCAAC
>11469752
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGCAGCGAGTTCCTCAC
"""
        self.assertEqual(obs, exp)

        with open(paths[2], 'U') as f:
            obs = f.read()
        exp = """>11468680
TAAACTGAAGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACG
"""
        self.assertEqual(obs, exp)

        with open(paths[3], 'U') as f:
            obs = f.read()
        exp = """>11458037
GACGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGGTTTCGAAGATCGGACTTCGAATTTCGAATTTCGAT
"""
        self.assertEqual(obs, exp)


class BufferedWriterTests(TestCase):

    def setUp(self):
        """ """
        self.files_to_remove = []
        tmp_dir = get_qiime_temp_dir()
        fd, self.test_fp = mkstemp(dir=tmp_dir,
                                   prefix='bufWriterTest',
                                   suffix='.txt')
        close(fd)
        self.files_to_remove.append(self.test_fp)

    def tearDown(self):
        """ """
        remove_files(self.files_to_remove)

    def test_init(self):
        """BufferedWriter constructor works"""

        BufferedWriter(self.test_fp)
        self.assertTrue(exists(self.test_fp))

    def test_write(self):
        """BufferedWriter writes nothing until max buffer reached."""

        b = BufferedWriter(self.test_fp, buf_size=2)
        b.write("1")
        content = open(self.test_fp, "r").readlines()
        self.assertEquals(content, [])

        # still nothing
        b.write("2")
        content = open(self.test_fp, "r").readlines()
        self.assertEquals(content, [])

        # finally, buffer is flushed
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

FILE1 = """>1 description field
ACCTACGTTAATACCCTGGTAGT
>2
ACCTACGTTAATACCCTGGTAGT
>3
AA"""

FILE2 = """>4 description field
ACCTACGTTAATACCCTGGTAGT
>5
ACCTACGTTAATACCCTGGTAGT
>6
AA
"""

FILE3 = """>7 description field
ACCTACGTTAATACCCTGGTAGT
>8
ACCTACGTTAATACCCTGGTAGT
>9
AA"""

LOGFILE1 = """Logging started at 12:17:52 on 01 Feb 2014
QIIME version: 1.8.0-dev, master@1d53d63"""

LOGFILE2 = """Logging started at 12:17:52 on 01 Feb 2014
QIIME version: 1.8.0-dev, master@1d53d63"""

LOGFILE3 = """Logging started at 12:17:52 on 01 Feb 2014
QIIME version: 1.8.0-dev, master@1d53d63"""

CONCAT_RESULT = """>1 description field
ACCTACGTTAATACCCTGGTAGT
>2
ACCTACGTTAATACCCTGGTAGT
>3
AA
>4 description field
ACCTACGTTAATACCCTGGTAGT
>5
ACCTACGTTAATACCCTGGTAGT
>6
AA
>7 description field
ACCTACGTTAATACCCTGGTAGT
>8
ACCTACGTTAATACCCTGGTAGT
>9
AA
"""

FASTA_FILE = """>11472286
GATGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGAACACTTTGTGTTTTGAGTTAATAGTTCGATAGTA
>11472384
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACGGGGGCAAC
>11468680
TAAACTGAAGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACG
>11458037
GACGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGGTTTCGAAGATCGGACTTCGAATTTCGAATTTCGAT
>11469739
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGAGAAGCTAACTTCTGA
>11469752
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGCAGCGAGTTCCTCAC
>11460523
AGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGCGAAATCGGGCACTCA
>11460543
TGGTTTGATCCTGGCTCAGGACAAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGAGAAGCCAGCTTTTGAT
>11480235
TGGTTTGATCCTGGCTCAGGATTAACGCTGGCGGCGCGCCTTATACATGCAAGTCGAACGAGCCTTGTGCTTCGCAC
>11480408
AATTTAGCGGCCGCGAATTCGCCCTTGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCA
"""

FASTA_FILE2 = """>11472286
GATGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGAACACTTTGTGTTTTGAGTTAATAGTTCGATAGTA
>11472384
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACGGGGGCAAC
>11468680
TAAACTGAAGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCACG
>11458037
GACGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGGTTTCGAAGATCGGACTTCGAATTTCGAATTTCGAT
>11469739
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGAGAAGCTAACTTCTGA
>11469752
AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGCAGCGAGTTCCTCAC
"""

if __name__ == "__main__":
    main()
