#!/usr/bin/env python

"""Tests of the trim_sff_primers.py file.

Note: this is presently just a stub that tests import and the one function that
isn't in the main block: this needs to be refactored.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", 'Kyle Bittinger']
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Release"

from cStringIO import StringIO
import os
import shutil
import tempfile

from cogent.util.unit_test import TestCase, main
from cogent.parse.binary_sff import (
    parse_binary_sff, write_binary_sff,
    )
from qiime.trim_sff_primers import (
    get_technical_lengths, set_sff_trimpoints, set_sff_trimpoints_with_sfftools,
    set_clip_qual_left, get_per_lib_sff_fps
    )


class TopLevelTests(TestCase):
    """Top-level tests of functions in trim_sff_primers"""

    def setUp(self):
        test_dir = os.path.dirname(os.path.abspath(__file__))
        self.sff_fn = 'F6AVWTA01.sff'
        orig_sff_fp = os.path.join(test_dir, 'test_support_files', 'F6AVWTA', self.sff_fn)
        self.sff_dir = tempfile.mkdtemp()
        self.sff_fp = os.path.join(self.sff_dir, self.sff_fn)
        shutil.copy(orig_sff_fp, self.sff_fp)

    def tearDown(self):
        shutil.rmtree(self.sff_dir)

    def test_get_technical_lengths(self):
        """get_technical_lengths should return correct dict of sample:length"""
        self.assertEqual(
            get_technical_lengths(mapping_with_linker), {'a':12, 'b':11})
        self.assertEqual(
            get_technical_lengths(mapping_without_linker), {'a':11, 'b':8})

    def test_set_sff_trimpoints_with_sfftools(self):
        _, orig_reads = parse_binary_sff(open(self.sff_fp), True)
        orig_reads = list(orig_reads)

        set_sff_trimpoints_with_sfftools(self.sff_dir, {'F6AVWTA01': 10})

        # check trimpoint file
        for line in open(self.sff_fp + '.trim'):
            toks = line.split()
            trim_start = int(toks[1])
            trim_end = int(toks[2])
            self.assertTrue(trim_start <= trim_end)
            self.assertEqual(trim_start, 11)

        # Check resultant SFF file
        _, reads = parse_binary_sff(open(self.sff_fp), True)
        for read, orig_read in zip(reads, orig_reads):
            self.assertEqual(read['clip_qual_left'], 11)
            # Check that eveything else is the same between original
            # reads and trimmed reads.
            orig_read['clip_qual_left'] = 11
            self.assertEqual(read, orig_read)

    def test_set_sff_trimpoints(self):
        _, orig_reads = parse_binary_sff(open(self.sff_fp), True)
        orig_reads = list(orig_reads)

        set_sff_trimpoints(self.sff_dir, {'F6AVWTA01': 10})

        _, reads = parse_binary_sff(open(self.sff_fp), True)
        for read, orig_read in zip(reads, orig_reads):
            self.assertEqual(read['clip_qual_left'], 11)
            # Check that eveything else is the same between original
            # reads and trimmed reads.
            orig_read['clip_qual_left'] = 11
            self.assertEqual(read, orig_read)

    def test_get_per_lib_sff_fps(self):
        self.assertEqual(
            list(get_per_lib_sff_fps(self.sff_dir)),
            [('F6AVWTA01', self.sff_fp)],
            )

    def test_set_clip_qual_left(self):
        orig_header, orig_reads = parse_binary_sff(open(self.sff_fp), True)
        orig_reads = list(orig_reads)

        _, clip_reads = set_clip_qual_left(
            (orig_header, orig_reads), 8)

        for read, orig_read in zip(clip_reads, orig_reads):
            self.assertEqual(read['clip_qual_left'], 9)
            # Check that eveything else is the same between original
            # reads and trimmed reads.
            orig_read['clip_qual_left'] = 9
            self.assertEqual(read, orig_read)


mapping_with_linker = StringIO("""\
#SampleID\tKEY_SEQ\tBARCODE\tLINKER\tPRIMER
a\tATGC\tCCC\tC\tCCCC
b\tATGC\tGG\tAAA\tCC
""")

mapping_without_linker = StringIO("""\
#SampleID\tKEY_SEQ\tBARCODE\tm\tPRIMER
a\tATGC\tCCC\tC\tCCCC
b\tATGC\tGG\tAAA\tCC
""")

if __name__ == '__main__':
    main()
