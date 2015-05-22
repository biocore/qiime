#!/usr/bin/env python
# file test_make_otu_table

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"  # consider project name
__credits__ = ["Rob Knight", "Justin Kuczynski", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import json
from unittest import TestCase, main
from StringIO import StringIO

from biom.table import Table
from biom.parse import parse_biom_table
import numpy as np

from qiime.make_otu_table import (libs_from_seqids, seqids_from_otu_to_seqid,
                                  make_otu_table)
from qiime.parse import parse_mapping_file, mapping_file_to_dict

class TopLevelTests(TestCase):

    """Tests of top-level functions"""

    def assertEqualOtuTable(self, obs, exp):
        """ """
        obs = json.loads(obs)
        exp = json.loads(exp)
        for e in ['generated_by', 'date']:
            del obs[e]
            del exp[e]
        self.assertEqual(obs, exp)

    def test_libs_from_seqids(self):
        """libs_from_seqids should identify correct libs"""
        seqids = ['ABC_001', 'DEF_002', 'ABC_003', 'GHI_JKL_001']
        self.assertEqual(libs_from_seqids(seqids),
                         set(['ABC', 'DEF', 'GHI_JKL']))

    def test_seqids_from_otu_to_seqid(self):
        """seqids_from_otu_to_seqid should return right seqids"""
        otu_to_seqid = {'0': ['ABC_0', 'DEF_1'], 'x': ['GHI_2']}
        self.assertEqual(seqids_from_otu_to_seqid(otu_to_seqid),
                         set(['ABC_0', 'DEF_1', 'GHI_2']))

    def test_make_otu_table_no_taxonomy(self):
        """make_otu_table should work without tax (new-style OTU table)"""
        otu_map_lines = """0	ABC_0	DEF_1
1	ABC_1
x	GHI_2	GHI_3	GHI_77
z	DEF_3	XYZ_1""".split('\n')

        obs = make_otu_table(otu_map_lines)
        data = [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 3, 0], [0, 1, 0, 1]]
        exp = Table(data, ['0', '1', 'x', 'z'], ['ABC', 'DEF', 'GHI', 'XYZ'],
                    input_is_dense=True)

        self.assertEqual(obs, exp)

    def test_make_otu_table_taxonomy(self):
        """make_otu_table should work with taxonomy"""
        otu_map_lines = """0	ABC_0	DEF_1
1	ABC_1
x	GHI_2	GHI_3	GHI_77
z	DEF_3	XYZ_1""".split('\n')
        taxonomy = {'0': ['Bacteria', 'Firmicutes'],
                    'x': ['Bacteria', 'Bacteroidetes']}
        obs = make_otu_table(otu_map_lines, taxonomy)

        data = [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 3, 0], [0, 1, 0, 1]]
        obs_md = [{'taxonomy': ['Bacteria', 'Firmicutes']},
                  {'taxonomy': ['None']},
                  {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                  {'taxonomy': ['None']}]
        exp = Table(data, ['0', '1', 'x', 'z'], ['ABC', 'DEF', 'GHI', 'XYZ'],
                    observation_metadata=obs_md, input_is_dense=True)

        self.assertEqual(obs, exp)

    def test_make_otu_table_with_sample_metadata(self):
        # Want to make sure that the order of the sample IDs in the OTU
        # map and the order of the IDs in the mapping file do not matter
        otu_map_lines = """0	ABC_0	DEF_1
1	ABC_1
x	GHI_2	GHI_3	GHI_77
z	DEF_3	XYZ_1""".split('\n')
        mapping_f = StringIO(MAPPING_FILE)
        sample_ids = ['ABC', 'DEF', 'GHI', 'XYZ']
        data = [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 3, 0], [0, 1, 0, 1]]

        map_data, map_header, map_comments = parse_mapping_file(mapping_f)
        sample_metadata = mapping_file_to_dict(map_data, map_header)

        sample_md = [sample_metadata[sample_id] for sample_id in sample_ids]

        obs = make_otu_table(otu_map_lines, sample_metadata=sample_metadata)
        exp = Table(data, ['0', '1', 'x', 'z'], sample_ids,
                    sample_metadata=sample_md, input_is_dense=True)

        self.assertEqual(obs, exp)

        # Test with a mapping file that is missing a sample's metadata,
        # make sure it raises the KeyError
        mapping_f = StringIO(MAPPING_FILE_MISSING_SAMPLE)
        map_data, map_header, map_comments = parse_mapping_file(mapping_f)
        sample_metadata = mapping_file_to_dict(map_data, map_header)

        with self.assertRaises(KeyError):
            obs = make_otu_table(otu_map_lines,
                                 sample_metadata=sample_metadata)


MAPPING_FILE = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Description
ABC	ATGC	AAAAAA	First Sample
XYZ	TGCA	AAAAAA	Fourth Sample
GHI	CATG	AAAAAA	Third Sample
DEF	GCAT	AAAAAA	Second Sample
"""

MAPPING_FILE_MISSING_SAMPLE = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Description
ABC	ATGC	AAAAAA	First Sample
XYZ	TGCA	AAAAAA	Fourth Sample
DEF	GCAT	AAAAAA	Second Sample
"""

if __name__ == '__main__':
    main()
