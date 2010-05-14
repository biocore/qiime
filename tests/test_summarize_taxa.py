#!/usr/bin/env python

"""Tests of code for adding taxa to OTU table"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Rob Knight"] 
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from qiime.summarize_taxa import make_new_summary_file, \
	add_summary_category_mapping

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_import(self):
        """empty test just verifies import of module"""
        pass
    def test_make_new_summary_file(self):
        """make_new_summary_file works
        """
        otu_table="""#Full OTU Counts
#OTU ID\ts1\ts2\ts3\ts4\tConsensus Lineage
0\t1\t0\t2\t4\tRoot;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
1\t1\t2\t0\t1\tRoot;Bacteria;Firmicutes;"Clostridia"
2\t0\t1\t1\t0\tRoot;Bacteria;Firmicutes;"Clostridia"
3\t1\t2\t1\t0\tRoot;Bacteria""".split('\n')
        result = make_new_summary_file(otu_table, 3, ';', 'False')
        self.assertEqual(result, ['#Full OTU Counts\n', 'Taxon\ts1\ts2\ts3\ts4\n', 'Root;Bacteria;Actinobacteria\t1.0\t0.0\t2.0\t4.0\n', 'Root;Bacteria;Firmicutes\t1.0\t3.0\t1.0\t1.0\n', 'Root;Bacteria;Other\t1.0\t2.0\t1.0\t0.0\n'])
        #test that works with relative abundances
        result = make_new_summary_file(otu_table, 3, ';', 'True')
        self.assertEqual(result, ['#Full OTU Counts\n', 'Taxon\ts1\ts2\ts3\ts4\n', 'Root;Bacteria;Actinobacteria\t0.333333333333\t0.0\t0.5\t0.8\n', 'Root;Bacteria;Firmicutes\t0.333333333333\t0.6\t0.25\t0.2\n', 'Root;Bacteria;Other\t0.333333333333\t0.4\t0.25\t0.0\n'])

    def test_add_summary_category_mapping(self):
        """make_new_summary_file works
        """
        otu_table="""#Full OTU Counts
#OTU ID\ts1\ts2\ts3\ts4\tConsensus Lineage
0\t1\t0\t2\t4\tRoot;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
1\t1\t2\t0\t1\tRoot;Bacteria;Firmicutes;"Clostridia"
2\t0\t1\t1\t0\tRoot;Bacteria;Firmicutes;"Clostridia"
3\t1\t2\t1\t0\tRoot;Bacteria""".split('\n')
        cat_mapping="""#SampleID\tBarcodeSequence\tTreatment\tDescription
#Test mapping file
s1\tAAAA\tControl\tControl mouse, I.D. 354
s2\tGGGG\tControl\tControl mouse, I.D. 355
s3\tCCCC\tExp\tDisease mouse, I.D. 356
s4\tTTTT\tExp\tDisease mouse, I.D. 357""".split('\n')
        result = add_summary_category_mapping(otu_table, cat_mapping, 3, ';',
            'False')
        self.assertEqual(result, ['#SampleID\tBarcodeSequence\tTreatment\tDescription\tRoot;Bacteria;Actinobacteria\tRoot;Bacteria;Firmicutes\tRoot;Bacteria;Other\n', '#Test mapping file\n', 's1\tAAAA\tControl\tControl mouse, I.D. 354\t1.0\t1.0\t1.0\n', 's2\tGGGG\tControl\tControl mouse, I.D. 355\t0.0\t3.0\t2.0\n', 's3\tCCCC\tExp\tDisease mouse, I.D. 356\t2.0\t1.0\t1.0\n', 's4\tTTTT\tExp\tDisease mouse, I.D. 357\t4.0\t1.0\t0.0\n'])
        result = add_summary_category_mapping(otu_table, cat_mapping, 3, ';',
            'True')
        self.assertEqual(result, ['#SampleID\tBarcodeSequence\tTreatment\tDescription\tRoot;Bacteria;Actinobacteria\tRoot;Bacteria;Firmicutes\tRoot;Bacteria;Other\n', '#Test mapping file\n', 's1\tAAAA\tControl\tControl mouse, I.D. 354\t0.333333333333\t0.333333333333\t0.333333333333\n', 's2\tGGGG\tControl\tControl mouse, I.D. 355\t0.0\t0.6\t0.4\n', 's3\tCCCC\tExp\tDisease mouse, I.D. 356\t0.5\t0.25\t0.25\n', 's4\tTTTT\tExp\tDisease mouse, I.D. 357\t0.8\t0.2\t0.0\n'])


#run unit tests if run from command-line
if __name__ == '__main__':
    main()
