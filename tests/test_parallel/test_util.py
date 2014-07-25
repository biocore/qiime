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
from tempfile import mkstemp
from unittest import TestCase, main

from skbio.util.misc import remove_files

from qiime.util import get_qiime_temp_dir
from qiime.parallel.util import (BufferedWriter, concatenate_files,
                                 merge_files_from_dirs, input_fasta_splitter)


class UtilTests(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_concatenate_files(self):
        pass

    def test_merge_files_from_dirs(self):
        pass

    def test_input_fasta_splitter(self):
        pass


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

if __name__ == "__main__":
    main()
