#!/usr/bin/env python
#file test_make_otu_table

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Rob Knight", "Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from qiime.make_otu_table import (libs_from_seqids,
        seqids_from_otu_to_seqid, make_otu_table, remove_otus)

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


    def test_make_otu_table_no_taxonomy_legacy(self):
        """make_otu_table should work without tax (legacy OTU table)"""
        otu_to_seqid ={ '0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        obs = make_otu_table(otu_to_seqid)
        exp = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ
0\t1\t1\t0\t0
1\t1\t0\t0\t0
x\t0\t0\t3\t0
z\t0\t1\t0\t1""" % __version__
        self.assertEqual(obs, exp)

    def test_make_otu_table_taxonomy_legacy(self):
        """make_otu_table should work wit tax (legacy OTU table)"""
        otu_to_seqid ={ '0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        taxonomy = {'0':'Bacteria;Firmicutes', 'x':'Bacteria;Bacteroidetes'}
        obs = make_otu_table(otu_to_seqid, taxonomy)
        exp = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes
1\t1\t0\t0\t0\tNone
x\t0\t0\t3\t0\tBacteria;Bacteroidetes
z\t0\t1\t0\t1\tNone""" % __version__
        self.assertEqual(obs, exp)

    def test_make_otu_table_no_taxonomy(self):
        """make_otu_table should work without tax (new-style OTU table)"""
        otu_to_seqid ={ '0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        obs = make_otu_table(otu_to_seqid, legacy=False)
        exp = """# QIIME v%s OTU table
OTU ID\tABC\tDEF\tGHI\tXYZ
0\t1\t1\t0\t0
1\t1\t0\t0\t0
x\t0\t0\t3\t0
z\t0\t1\t0\t1""" % __version__
        self.assertEqual(obs, exp)

    def test_make_otu_table_taxonomy(self):
        """make_otu_table should work wit tax (new-style OTU table)"""
        otu_to_seqid ={ '0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        taxonomy = {'0':'Bacteria;Firmicutes', 'x':'Bacteria;Bacteroidetes'}
        obs = make_otu_table(otu_to_seqid, taxonomy, legacy=False)
        exp = """# QIIME v%s OTU table
OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes
1\t1\t0\t0\t0\tNone
x\t0\t0\t3\t0\tBacteria;Bacteroidetes
z\t0\t1\t0\t1\tNone""" % __version__
        self.assertEqual(obs, exp)
        
    def test_remove_otus(self):
        """remove_otus functions as expected """
        otu_to_seqid = {'0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        otus_to_exclude = ['0 some comment','42 not a real otu id','z']
        expected = {'1':['ABC_1'],
                    'x':['GHI_2', 'GHI_3','GHI_77']}
        actual = remove_otus(otu_to_seqid,otus_to_exclude)
        self.assertEqual(actual,expected)
        
        # functions with empty otus_to_exclude
        otu_to_seqid = {'0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        expected = {'0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        actual = remove_otus(otu_to_seqid,[])
        self.assertEqual(actual,expected)
        
        # functions with complete otus_to_exclude
        otu_to_seqid = {'0':['ABC_0','DEF_1'],
                        '1':['ABC_1'],
                        'x':['GHI_2', 'GHI_3','GHI_77'],
                        'z':['DEF_3','XYZ_1']
                        }
        actual = remove_otus(otu_to_seqid,list('01xz'))
        self.assertEqual(actual,{})
        
        


if __name__ =='__main__':
    main()
