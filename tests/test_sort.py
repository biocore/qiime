#!/usr/bin/env python
# File created on 15 Feb 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.parse import parse_otu_table, parse_mapping_file
from qiime.sort import (sort_sample_ids_by_mapping_value,
                        sort_fasta_by_abundance, natsort, sort_otu_table,
                        sort_otu_table_by_mapping_field)

class SortTests(TestCase):
    
    def setUp(self):
        self.mapping_f1 = mapping_f1.split('\n')
        self.otu_table1 = otu_table1.split('\n')
        self.mapping_f2 = mapping_f2.split('\n')
        self.age_sorted_otu_table1 = age_sorted_otu_table1.split('\n')
        self.name_sorted_otu_table1 = name_sorted_otu_table1.split('\n')
        self.nothing_sorted_otu_table1 = nothing_sorted_otu_table1.split('\n')
        self.otu_table1_bad_sampleID = otu_table1_bad_sampleID.split('\n')
        self.dirs_to_remove = []
        self.files_to_remove = []

    def tearDown(self):
        for dir in  self.dirs_to_remove:
            if exists(dir):
                rmdir(dir)
        remove_files(self.files_to_remove)
    
    def test_sort_sample_ids_by_mapping_value(self):
        """ sort_sample_ids_by_mapping_value functions as expected """
        actual = sort_sample_ids_by_mapping_value(mapping_file=self.mapping_f1,
                                         field='days_since_epoch',
                                         field_type_f=float)
        expected = zip(['NotInOtuTable','1','Z2','Z1','A'],
                       [0.0,5.7,10,23,400000])
        self.assertEqual(actual,expected)
        
    def test_sort_sample_ids_by_mapping_value_error(self):
        """ sort_sample_ids_by_mapping_value handles errors """
        self.assertRaises(ValueError,
                          sort_sample_ids_by_mapping_value,
                          mapping_file=self.mapping_f1,
                          field='years_since_spoch',
                          field_type_f=float)
                          
        self.assertRaises(ValueError,
                          sort_sample_ids_by_mapping_value,
                          mapping_file=self.mapping_f1,
                          field='Something',
                          field_type_f=float)
                          
    def test_sort_fasta_by_abundance(self):
        """sort_fasta_by_abundance functions as expected"""
        class FakeOutF(object):
          def __init__(self):
              self.s = ""
          def write(self,s):
              self.s += s

        actual = FakeOutF()
        expected = ""
        sort_fasta_by_abundance([],actual)
        self.assertEqual(actual.s,expected)

        # no sorting necessary
        actual = FakeOutF()
        expected1 = "\n".join(['>s1','ACCGT','>s2 comment','ATTGC',''])
        expected2 = "\n".join(['>s2 comment','ATTGC','>s1','ACCGT',''])
        sort_fasta_by_abundance(['>s1','ACCGT','>s2 comment','ATTGC'],actual)
        # order is unimportant here
        self.assertTrue(actual.s == expected1 or actual.s == expected2)

        # sorting necessary
        actual = FakeOutF()
        inseqs = ['>s1','ACCGT',
                 '>s2 comment','ATTGC',
                 '>s3 blah','ATTGC']
        expected = "\n".join(['>s2 comment','ATTGC',
                            '>s3 blah','ATTGC',
                            '>s1','ACCGT',''])
        sort_fasta_by_abundance(inseqs,actual)
        self.assertEqual(actual.s,expected)


    def test_natsort(self):
        """natsort should perform numeric comparisons on strings"""
        # string with alpha and numerics sort correctly
        s = 'sample1 sample2 sample11 sample12'.split()
        self.assertEqual(natsort(s), 
          'sample1 sample2 sample11 sample12'.split())
        s.reverse()
        self.assertEqual(natsort(s), 
          'sample1 sample2 sample11 sample12'.split())
        self.assertEqual(natsort(list('cba321')),list('123abc'))

        # strings with alpha only sort correctly
        self.assertEqual(natsort(list('cdba')),list('abcd'))

        # string of ints sort correctly
        self.assertEqual(natsort(['11','2','1','0']),
                               ['0','1','2','11'])

        # strings of floats sort correctly
        self.assertEqual(natsort(['1.11','1.12','1.00','0.009']),
                               ['0.009','1.00','1.11','1.12'])

        # string of ints sort correctly
        self.assertEqual(natsort([('11','A'),('2','B'),('1','C'),('0','D')]),
                            [('0','D'),('1','C'),('2','B'),('11','A')])

    def test_sort_otu_table_by_mapping_field_all_values_differ(self):
        """ sort_otu_table_by_mapping_field fns when all values differ"""

        actual = sort_otu_table_by_mapping_field(parse_otu_table(self.otu_table1),
                                parse_mapping_file(self.mapping_f2),
                                sort_field = "Age")
        expected = parse_otu_table(self.age_sorted_otu_table1)
        # sample ids match expected
        self.assertEqual(actual[0],expected[0])
        # otu ids match expected
        self.assertEqual(actual[1],expected[1])
        # otu data match expected
        self.assertEqual(actual[2],expected[2])
        # taxa match expected
        self.assertEqual(actual[3],expected[3])
        
    def test_sort_otu_table(self):
        """ sort_otu_table fns as expected """

        actual = sort_otu_table(parse_otu_table(self.otu_table1),
                                ['NA','Key','Fing'])
        expected = parse_otu_table(self.age_sorted_otu_table1)
        # sample ids match expected
        self.assertEqual(actual[0],expected[0])
        # otu ids match expected
        self.assertEqual(actual[1],expected[1])
        # otu data match expected
        self.assertEqual(actual[2],expected[2])
        # taxa match expected
        self.assertEqual(actual[3],expected[3])

    def test_sort_otu_table_error(self):
        """ sort_otu_table handles errors """

        self.assertRaises(ValueError,sort_otu_table,
            parse_otu_table(self.otu_table1),['NA','Key','Fing','Key'])
        self.assertRaises(KeyError,sort_otu_table,
            parse_otu_table(self.otu_table1),['NA','Key'])

    def test_sort_otu_table_by_mapping_field_some_values_differ(self):
        """ sort_otu_table fns when some values differ"""

        actual = sort_otu_table_by_mapping_field(parse_otu_table(self.otu_table1),
                              parse_mapping_file(self.mapping_f2),
                              sort_field = "Nothing")
        expected = parse_otu_table(self.nothing_sorted_otu_table1)
        # sample ids match expected
        self.assertEqual(actual[0],expected[0])
        # otu ids match expected
        self.assertEqual(actual[1],expected[1])
        # otu data match expected
        self.assertEqual(actual[2],expected[2])
        # taxa match expected
        self.assertEqual(actual[3],expected[3])

    def test_sort_otu_table_by_mapping_field_some_values_same(self):
        """ sort_otu_table_by_mapping_field fns when all values are the same"""

        actual = sort_otu_table_by_mapping_field(parse_otu_table(self.otu_table1),
                              parse_mapping_file(self.mapping_f2),
                              sort_field = "Name")
        expected = parse_otu_table(self.name_sorted_otu_table1)
        # sample ids match expected
        self.assertEqual(actual[0],expected[0])
        # otu ids match expected
        self.assertEqual(actual[1],expected[1])
        # otu data match expected
        self.assertEqual(actual[2],expected[2])
        # taxa match expected
        self.assertEqual(actual[3],expected[3])

    def test_sort_otu_table_by_mapping_field_error(self):
        """ sort_otu_table_by_mapping_field fails on samples in otu table but not mapping"""

        self.assertRaises(KeyError,sort_otu_table_by_mapping_field,
                                   parse_otu_table(self.otu_table1_bad_sampleID),
                                   parse_mapping_file(self.mapping_f2),
                                   sort_field = "Age")



mapping_f1 = """#SampleID\tSomething\tdays_since_epoch
Z1\t42\t23
Z2\thello\t10
A\t4\t400000
1\tr\t5.7
NotInOtuTable\tf\t0"""

otu_table1 = """# Some comment
OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
# comment
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

otu_table1_bad_sampleID = """# Some comment
OTU ID	Fing	Key	NotInMapping	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
# comment
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

age_sorted_otu_table1 = """# Some comment
OTU ID	NA	Key	Fing	Consensus Lineage
0	42	44536	19111	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	6	3500	1216	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	2	1184	1803	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	17	4903	1722	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	34	2074	589	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

nothing_sorted_otu_table1 = """# Some comment
OTU ID	Fing	NA	Key	Consensus Lineage
0	19111	42	44536	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	6	3500	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	2	1184	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	17	4903	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	34	2074	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

name_sorted_otu_table1 = """# Some comment
OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

# values in 'Age' column sort differently incorrectly if sorted as strings
mapping_f2 = """
#SampleID	Name	Age	Nothing
Fing	Blah	11	4
Key	Blah	2	5
NA	Blah	1	4
IgnoredSample	Blah	1	4"""


if __name__ == "__main__":
    main()