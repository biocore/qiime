#!/usr/bin/env python
#file test_filter_by_metadata.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

from qiime.parse import parse_otu_table, parse_mapping_file, parse_metadata_state_descriptions
from sys import argv
from string import strip
from cogent.util.unit_test import TestCase, main
from numpy import array
from StringIO import StringIO
from qiime.filter_by_metadata import (get_sample_ids, find_good_cols,
    filter_line, filter_map, filter_otus_and_map)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """Define a few simple tables"""
        self.otu_str = """#Full OTU Counts
#OTU ID\ta\tb\tc\td\te
1\t1\t2\t3\t4\t5
2\t5\t4\t3\t2\t1"""
        self.otu_table = parse_otu_table(StringIO(self.otu_str))
        self.otu_tax_str = """#Full OTU Counts
#OTU ID\ta\tb\tc\td\te\tConsensus Lineage
1\t1\t2\t3\t4\t5\tBacteria:Firmicutes
2\t5\t4\t3\t2\t1\tBacteria:Proteobacteria"""
        self.otu_tax_table = parse_otu_table(StringIO(self.otu_tax_str))
        self.map_str = """#SampleID\tStudy\tBodySite\tDescription
a\tDog\tStool\tx
b\tDog\tStool\ty
c\tHand\tPalm\tz
d\tWholeBody\tPalm\ta
e\tWholeBody\tStool\tb"""
        self.map_data, self.map_headers, self.map_comments =\
         parse_mapping_file(StringIO(self.map_str))

    def test_get_sample_ids(self):
        """get_sample_ids should return sample ids matching criteria."""
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,\
            parse_metadata_state_descriptions('Study:Twin')), [])
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,\
            parse_metadata_state_descriptions('Study:Dog')), ['a','b'])
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,\
            parse_metadata_state_descriptions('Study:*,!Dog')), ['c','d','e'])
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,\
            parse_metadata_state_descriptions('Study:*,!Dog;BodySite:Stool')), ['e'])
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,\
            parse_metadata_state_descriptions('BodySite:Stool')), ['a','b','e'])

    def test_find_good_cols(self):
        """find_good_cols should return col indices of specified samples."""
        map_line = self.otu_str.splitlines()[1]
        map_tax_line = self.otu_tax_str.splitlines()[1]
        self.assertEqual(find_good_cols(map_line, ['b','e']), [0,2,5])
        #note: always returns col 0 for the OTU id
        self.assertEqual(find_good_cols(map_tax_line, ['b','e']), [0,2,5,-1])
        #note: adds -1 for the taxonomy col

    def test_filter_line(self):
        """filter_line should return correct trimmed OTU line."""
        result = StringIO('')
        line = self.otu_tax_str.splitlines()[-1]
        filter_line(line, [0,2,5,-1], None, result)
        result.seek(0)
        self.assertEqual(result.read(), '2\t4\t1\tBacteria:Proteobacteria\n')
        #test without taxonomy
        result = StringIO('')
        line = self.otu_str.splitlines()[-1]
        filter_line(line, [0,2,5], None, result)
        result.seek(0)
        self.assertEqual(result.read(), '2\t4\t1\n')
        #test with min count
        result = StringIO('')
        line = self.otu_str.splitlines()[-1]
        filter_line(line, [0,1,2,3], 20, result)
        result.seek(0)
        self.assertEqual(result.read(), '')
        filter_line(line, [0,1,2,3], 5, result)
        result.seek(0)
        self.assertEqual(result.read(), '2\t5\t4\t3\n')

    def test_filter_map(self):
        """filter_map should filter map file according to sample ids"""
        self.assertEqual(filter_map(self.map_data, self.map_headers,\
         ['a','b','c','d','e']), (self.map_headers, self.map_data))
        self.assertEqual(filter_map(self.map_data, self.map_headers, ['a']),
            (['SampleID','Description'],['a\tx'.split('\t')]))

    def test_filter_otus_and_map(self):
        """filter_otus_and_map should filter both OTU and map files."""
        def get_result(valid_states_str, num_seqs_per_otu):
            otu_infile = StringIO(self.otu_tax_str)
            map_infile = StringIO(self.map_str)
            otu_outfile = StringIO('')
            map_outfile = StringIO('')
            filter_otus_and_map(map_infile, otu_infile, map_outfile, otu_outfile,
                valid_states_str, num_seqs_per_otu)
            otu_outfile.seek(0)
            map_outfile.seek(0)
            return otu_outfile.read(), map_outfile.read()
        self.assertEqual(get_result('Study:*', None), 
            (self.otu_tax_str+'\n', self.map_str))

        no_dog_otu="""#Full OTU Counts
#OTU ID\tc\td\te\tConsensus Lineage
1\t3\t4\t5\tBacteria:Firmicutes
2\t3\t2\t1\tBacteria:Proteobacteria"""
        no_dog_map = """#SampleID\tStudy\tBodySite\tDescription
c\tHand\tPalm\tz
d\tWholeBody\tPalm\ta
e\tWholeBody\tStool\tb"""

        self.assertEqual(get_result('Study:*,!Dog', None),
            (no_dog_otu+'\n', no_dog_map))




#run tests if called from command line
if __name__ == "__main__":
    main()
