#!/usr/bin/env python
#file test_make_otu_table

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Rob Knight", "Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from qiime.make_otu_table import (libs_from_seqids,
        seqids_from_otu_to_seqid, make_otu_map)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""
    def test_libs_from_seqids(self):
        """libs_from_seqids should identify correct libs"""
        seqids = ['ABC_001', 'DEF_002', 'ABC_003', 'GHI_JKL_001']
        self.assertEqual(libs_from_seqids(seqids),
                set(['ABC', 'DEF', 'GHI_JKL']))

    def test_seqids_from_otu_to_seqid(self):
        """seqids_from_otu_to_seqid should return right seqids"""
        otu_to_seqid ={'0':['ABC_0','DEF_1'],'x':['GHI_2']}
        self.assertEqual(seqids_from_otu_to_seqid(otu_to_seqid),
            set(['ABC_0', 'DEF_1', 'GHI_2']))

    def test_make_otu_map_no_taxonomy(self):
        """make_otu_map should work without supplied taxonomy"""
        otu_to_seqid ={ '0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        obs = make_otu_map(otu_to_seqid)
        exp = """#Full OTU Counts
#OTU ID\tABC\tDEF\tGHI\tXYZ
0\t1\t1\t0\t0
1\t1\t0\t0\t0
x\t0\t0\t3\t0
z\t0\t1\t0\t1"""
        self.assertEqual(obs, exp)

    def test_make_otu_map_taxonomy(self):
        """make_otu_map should work with supplied taxonomy"""
        otu_to_seqid ={ '0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        taxonomy = {'0':'Bacteria;Firmicutes', 'x':'Bacteria;Bacteroidetes'}
        obs = make_otu_map(otu_to_seqid, taxonomy)
        exp = """#Full OTU Counts
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes
1\t1\t0\t0\t0\tNone
x\t0\t0\t3\t0\tBacteria;Bacteroidetes
z\t0\t1\t0\t1\tNone"""
        self.assertEqual(obs, exp)


if __name__ =='__main__':
    main()
