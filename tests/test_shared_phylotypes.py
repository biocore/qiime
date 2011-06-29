#!/usr/bin/env python

"""Tests for computing shared phylotypes."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_otu_table
from qiime.shared_phylotypes import _calc_shared_phylotypes_pairwise,\
    _calc_shared_phylotypes_multiple, calc_shared_phylotypes

class Test_shared_phylotypes(TestCase):
    def setUp(self):
        self.otu_table_as_string = ["#Test otu table",
                       "\t".join(["#OTU ID","S1","S2","S3"]),
                       "\t".join(["0",      "1" ,"0" ,"2" ]),
                       "\t".join(["1",      "1" ,"2" ,"0" ]),
                       "\t".join(["2",      "1" ,"0" ,"0" ]),
                       "\t".join(["3",      "1" ,"0" ,"2" ]),
                       "\t".join(["4",      "1" ,"1" ,"2" ])]

        _, _, self.otu_table, _ = parse_otu_table(self.otu_table_as_string)

    def test_calc_shared_phylotypes_pairwise(self):
        """_calc_shared_phylotypes_pairwise works as expected"""

        self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 0, 0), 5)
        self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 0, 1), 2)
        self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 0, 2), 3)
        self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 2, 2), 3)

    def test_calc_shared_phylotypes_multiple(self):
        """_calc_shared_phylotypes_multiple works as expected"""

        #test for <2 idxes
        self.assertRaises(ValueError, _calc_shared_phylotypes_multiple, self.otu_table, [])
        self.assertRaises(ValueError, _calc_shared_phylotypes_multiple, self.otu_table, [0])

        #test that func is identical to _calc_shared_phylotypes_pairwise with 2 idx
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, [0, 0]), 5)
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, [0, 1]), 2)
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, [0, 2]), 3)
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, [2, 2]), 3)

        # works with more than 2 samples
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, [0,1,2]), 1)


    def test_calc_shared_phylotypes(self):
        """calc_shared_phylotypes computes correct matrix"""
        
        observed = calc_shared_phylotypes(iter(self.otu_table_as_string))
        expected = """\tS1\tS2\tS3
S1\t5\t2\t3
S2\t2\t2\t1
S3\t3\t1\t3\n"""
        self.assertEqual(observed, expected)
                       
                                 
if __name__ == "__main__":
    main()
