#!/usr/bin/env python

"""Tests for computing shared phylotypes."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder","Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from biom.parse import parse_biom_table_str
from qiime.shared_phylotypes import _calc_shared_phylotypes_pairwise,\
    _calc_shared_phylotypes_multiple, calc_shared_phylotypes

class Test_shared_phylotypes(TestCase):
    def setUp(self):
        #self.otu_table_as_string = ["#Test otu table",
        #               "\t".join(["#OTU ID","S1","S2","S3"]),
        #               "\t".join(["0",      "1" ,"0" ,"2" ]),
        #               "\t".join(["1",      "1" ,"2" ,"0" ]),
        #               "\t".join(["2",      "1" ,"0" ,"0" ]),
        #               "\t".join(["3",      "1" ,"0" ,"2" ]),
        #               "\t".join(["4",      "1" ,"1" ,"2" ])]
        self.biom_as_string = '{"rows": [{"id": "0", "metadata": null}, {"id": "1", "metadata": null}, {"id": "2", "metadata": null}, {"id": "3", "metadata": null}, {"id": "4", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 2, 2.0], [1, 0, 1.0], [1, 1, 2.0], [2, 0, 1.0], [3, 0, 1.0], [3, 2, 2.0], [4, 0, 1.0], [4, 1, 1.0], [4, 2, 2.0]], "columns": [{"id": "S1", "metadata": null}, {"id": "S2", "metadata": null}, {"id": "S3", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2583", "matrix_type": "sparse", "shape": [5, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-22T01:06:31.645277", "type": "OTU table", "id": null, "matrix_element_type": "float"}'
        self.otu_table = parse_biom_table_str(self.biom_as_string)

    def test_calc_shared_phylotypes_pairwise(self):
        """_calc_shared_phylotypes_pairwise works as expected"""

        #self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 0, 0), 5)
        #self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 0, 1), 2)
        #self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 0, 2), 3)
        #self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 2, 2), 3)
        self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 'S1', 'S1'), 5)
        self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 'S1', 'S2'), 2)
        self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 'S1', 'S3'), 3)
        self.assertEqual(_calc_shared_phylotypes_pairwise(self.otu_table, 'S3', 'S3'), 3)

    def test_calc_shared_phylotypes_multiple(self):
        """_calc_shared_phylotypes_multiple works as expected"""

        #test for <2 idxes
        self.assertRaises(ValueError, _calc_shared_phylotypes_multiple, self.otu_table, [])
        self.assertRaises(ValueError, _calc_shared_phylotypes_multiple, self.otu_table, ['S1'])

        #test that func is identical to _calc_shared_phylotypes_pairwise with 2 idx
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, ['S1', 'S1']), 5)
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, ['S1', 'S2']), 2)
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, ['S1', 'S3']), 3)
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, ['S3', 'S3']), 3)

        # works with more than 2 samples
        self.assertEqual(_calc_shared_phylotypes_multiple(self.otu_table, ['S1','S2','S3']), 1)


    def test_calc_shared_phylotypes(self):
        """calc_shared_phylotypes computes correct matrix"""
        
        observed = calc_shared_phylotypes(self.biom_as_string)
        expected = """\tS1\tS2\tS3
S1\t5\t2\t3
S2\t2\t2\t1
S3\t3\t1\t3\n"""
        self.assertEqual(observed, expected)
                       
                                 
if __name__ == "__main__":
    main()
