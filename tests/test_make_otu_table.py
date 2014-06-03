#!/usr/bin/env python
# file test_make_otu_table

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"  # consider project name
__credits__ = ["Rob Knight", "Justin Kuczynski"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import json
from unittest import TestCase, main
from qiime.make_otu_table import (libs_from_seqids,
                                  seqids_from_otu_to_seqid, make_otu_table)
from biom.table import Table
from biom.parse import parse_biom_table


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

        obs = make_otu_table(otu_map_lines, constructor=Table)
        exp = """{"rows": [{"id": "0", "metadata": null}, {"id": "1", "metadata": null}, {"id": "x", "metadata": null}, {"id": "z", "metadata": null}], "format": "Biological Observation Matrix 0.9dev", "data": [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 3, 0], [0, 1, 0, 1]], "columns": [{"id": "ABC", "metadata": null}, {"id": "DEF", "metadata": null}, {"id": "GHI", "metadata": null}, {"id": "XYZ", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2532", "matrix_type": "dense", "shape": [4, 4], "format_url": "http://biom-format.org", "date": "2011-12-21T00:49:15.978315", "type": "OTU table", "id": null, "matrix_element_type": "int"}"""

        self.assertEqual(
            parse_biom_table(obs.split('\n'), input_is_dense=False),
            parse_biom_table(exp.split('\n'), input_is_dense=True))

    def test_make_otu_table_taxonomy(self):
        """make_otu_table should work with taxonomy"""
        otu_map_lines = """0	ABC_0	DEF_1
1	ABC_1
x	GHI_2	GHI_3	GHI_77
z	DEF_3	XYZ_1""".split('\n')
        taxonomy = {'0': ['Bacteria', 'Firmicutes'],
                    'x': ['Bacteria', 'Bacteroidetes']}
        obs = make_otu_table(
            otu_map_lines,
            taxonomy,
            constructor=Table)
        exp = """{"rows": [{"id": "0", "metadata": {"taxonomy": ["Bacteria", "Firmicutes"]}}, {"id": "1", "metadata": {"taxonomy": ["None"]}}, {"id": "x", "metadata": {"taxonomy": ["Bacteria", "Bacteroidetes"]}}, {"id": "z", "metadata": {"taxonomy": ["None"]}}], "format": "Biological Observation Matrix 0.9dev", "data": [[1.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 3.0, 0.0], [0.0, 1.0, 0.0, 1.0]], "columns": [{"id": "ABC", "metadata": null}, {"id": "DEF", "metadata": null}, {"id": "GHI", "metadata": null}, {"id": "XYZ", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2532", "matrix_type": "dense", "shape": [4, 4], "format_url": "http://biom-format.org", "date": "2011-12-21T00:19:30.961477", "type": "OTU table", "id": null, "matrix_element_type": "int"}"""

        self.assertEqual(
            parse_biom_table(obs.split('\n'), input_is_dense=False),
            parse_biom_table(exp.split('\n'), input_is_dense=True))


if __name__ == '__main__':
    main()
