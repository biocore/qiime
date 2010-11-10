#!/usr/bin/env python

from cStringIO import StringIO
import os
import shutil
import tempfile

from cogent.util.unit_test import TestCase, main
from qiime.sra_spreadsheet_to_map_files import (
    strip_quotes, collect_study_groups, remap_lines, write_map_files,
    get_study_groups,
    )

"""Tests of the sra_spreadsheet_to_map_files.py file.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Development"

class TopLevelTests(TestCase):
    """Top-level tests of functions in sra_spreadsheet_to_map_files.py"""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_strip_quotes(self):
        """strip_quotes should strip terminal quotes from field."""
        self.assertEqual(strip_quotes('abc'), 'abc')
        self.assertEqual(strip_quotes('"abc"'), 'abc')

    def test_remap_lines(self):
        """remap_lines should fix some issues with input lines."""
        lines = """POOL_MEMBER_NAME\tBARCODE\tPRIMER\tLINKER\tabc
x\tAA\tGGG\tCC\tx_x
y\tAC\tCCC\tAA\ty_y""".splitlines()
        result = remap_lines(lines[0].split('\t'), 
            [i.split('\t') for i in lines[1:]])
        self.assertEqual(result,
        [['#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence'] + 
            'POOL_MEMBER_NAME\tBARCODE\tPRIMER\tLINKER\tabc'.split('\t') +
            ['Description'],
            ['x','AA','CCGGG','x','AA','GGG','CC','x_x','None'],
            ['y','AC','AACCC','y','AC','CCC','AA','y_y','None']])

    def test_get_study_groups(self):
        obs_header, obs_groups = get_study_groups(StringIO("#" + experiment_txt))
        self.assertEqual(obs_header, [
            "EXPERIMENT_TITLE", "EXPERIMENT_CENTER", "STUDY_REF",
            "SAMPLE_ALIAS", "BARCODE", "PRIMER", "LINKER",
            "RUN_PREFIX", "PRIMER_READ_GROUP_TAG",
            ])
        self.assertEqual(obs_groups[('study1', 'REGION01')], [[
            'experiment1', 'center1', 'study1', 'sample1', 'AA', 'CC', 'GGG',
            'REGION01', 'A',
            ]])

    def test_write_map_files(self):
        input_fp = os.path.join(self.temp_dir, 'experiment.txt')
        open(input_fp, 'w').write(experiment_txt)
        write_map_files(input_fp)
        observed = os.listdir(self.temp_dir)
        observed.sort()
        expected = [
            'experiment.txt', 'study1_REGION01.map', 'study1_REGION02.map']
        self.assertEqual(observed, expected)

        observed_map = open(os.path.join(self.temp_dir, observed[1])).read()
        self.assertEqual(observed_map, map_txt)
 
experiment_txt = """\
EXPERIMENT_TITLE	EXPERIMENT_CENTER	STUDY_REF	SAMPLE_ALIAS	BARCODE	PRIMER	LINKER	RUN_PREFIX	PRIMER_READ_GROUP_TAG
experiment1	center1	study1	sample1	AA	CC	GGG	REGION01	A
experiment1	center1	study1	sample1	AA	CC	GGG	REGION02	A
"""

map_txt = """\
#SampleID	BarcodeSequence	LinkerPrimerSequence	EXPERIMENT_TITLE	EXPERIMENT_CENTER	STUDY_REF	SAMPLE_ALIAS	BARCODE	PRIMER	LINKER	RUN_PREFIX	PRIMER_READ_GROUP_TAG	REGION	EXPERIMENT_ALIAS	RUN_ALIAS	BARCODE_READ_GROUP_TAG	POOL_MEMBER_NAME	POOL_MEMBER_FILENAME	RUN_CENTER	STUDY_CENTER	SAMPLE_CENTER	DEFAULT_SAMPLE_CENTER	DEFAULT_SAMPLE_NAME	DEFAULT_SAMPLE_FILENAME	DEFAULT_RUN_ALIAS	PLATFORM	KEY_SEQ	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	Description
REGION01_sample1_A	AA	GGGCC	experiment1	center1	study1	sample1	AA	CC	GGG	REGION01	A	0	study1_REGION01	study1_sample1_REGION01	REGION01_AA	REGION01_sample1_A	REGION01_sample1_A.sff	center1	center1	center1	center1	study1_default	study1_default_REGION01.sff	study1_default_REGION01	Titanium	TCAG	AMPLICON	GENOMIC	PCR	None
"""

if __name__ == '__main__':
    main()
