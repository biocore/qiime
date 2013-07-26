#!/usr/bin/env python

"""Tests of code for adding taxa to OTU table"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Daniel McDonald", "Antonio Gonzalez Pena", "Jose Carlos Clemente Litran"] 
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from qiime.summarize_taxa import make_summary, \
	add_summary_mapping, sum_counts_by_consensus
from qiime.parse import parse_mapping_file
from qiime.util import convert_otu_table_relative
from numpy import array
from biom.table import table_factory
from biom.parse import parse_biom_table

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        self.otu_table_vals = array([[1,0,2,4],
                               [1,2,0,1],
                               [0,1,1,0],
                               [1,2,1,0]])
        
        {(0, 0):1.0, (0, 2):2.0, (0, 3):4.0,
                               (1, 0):1.0, (1, 1):2.0, (1, 3):1.0,
                               (2, 1):1.0, (2, 2):1.0, (3, 0):1.0,
                               (3, 1): 2.0, (3, 2):1.0}

        self.otu_table = table_factory(self.otu_table_vals,
                                        ['s1', 's2', 's3', 's4'],
                                        ['0', '1', '2', '3'],
                                        None,
                                        [{"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]},
                                         {"taxonomy": ["Root", "Bacteria", "Firmicutes", "\"Clostridia\""]},
                                         {"taxonomy": ["Root", "Bacteria", "Firmicutes", "\"Clostridia\""]},
                                         {"taxonomy": ["Root", "Bacteria"]}])

#        self.otu_table="""#Full OTU Counts
##OTU ID\ts1\ts2\ts3\ts4\tConsensus Lineage
#0\t1\t0\t2\t4\tRoot;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
#1\t1\t2\t0\t1\tRoot;Bacteria;Firmicutes;"Clostridia"
#2\t0\t1\t1\t0\tRoot;Bacteria;Firmicutes;"Clostridia"
#3\t1\t2\t1\t0\tRoot;Bacteria""".split('\n')

        self.mapping="""#SampleID\tBarcodeSequence\tTreatment\tDescription
#Test mapping file
s1\tAAAA\tControl\tControl mouse, I.D. 354
s2\tGGGG\tControl\tControl mouse, I.D. 355
s3\tCCCC\tExp\tDisease mouse, I.D. 356
s4\tTTTT\tExp\tDisease mouse, I.D. 357""".split('\n')

    def test_sum_counts_by_consensus(self):
        """should sum otu counts by consensus"""
        #otu_table = parse_otu_table(self.otu_table)
        #otu_table = parse_biom_table(self.otu_table)
        obs_result, obs_mapping = sum_counts_by_consensus(self.otu_table, 3)
        exp_result = {('Root','Bacteria','Actinobacteria'):array([1,0,2,4]),
                      ('Root','Bacteria','Firmicutes'):array([1,3,1,1]),
                      ('Root','Bacteria','Other'):array([1,2,1,0])}
        exp_mapping = {'s1':0, 's2':1, 's3':2, 's4':3}
        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obs_mapping, exp_mapping)

        obs_result, obs_mapping = sum_counts_by_consensus(self.otu_table, 2)
        exp_result = {('Root','Bacteria'):array([3,5,4,5])}
        exp_mapping = {'s1':0, 's2':1, 's3':2, 's4':3}
        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obs_mapping, exp_mapping)

        obs_result, obs_mapping = sum_counts_by_consensus(self.otu_table, 4)
        exp_result = {('Root','Bacteria','Actinobacteria','Actinobacteria'):\
                array([1,0,2,4]),
                      ('Root','Bacteria','Firmicutes','"Clostridia"'):\
                              array([1,3,1,1]),
                      ('Root','Bacteria','Other','Other'):array([1,2,1,0])}
        exp_mapping = {'s1':0, 's2':1, 's3':2, 's4':3}
        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obs_mapping, exp_mapping)

    def test_make_new_summary_file(self):
        """make_new_summary_file works
        """
        lower_percentage, upper_percentage = None, None
        #otu_table = parse_otu_table(self.otu_table, int)
        #otu_table = parse_biom_table(self.otu_table)
        summary, header = make_summary(self.otu_table, 3, upper_percentage, lower_percentage)
        self.assertEqual(header, ['Taxon','s1','s2','s3','s4'])
        self.assertEqual(summary, [[('Root','Bacteria','Actinobacteria'),1,0,2,4], 
                                   [('Root','Bacteria','Firmicutes'),1,3,1,1], 
                                   [('Root','Bacteria','Other'),1,2,1,0]])

        #test that works with relative abundances
        #otu_table = parse_otu_table(self.otu_table, float)
        #otu_table = parse_biom_table(self.otu_table, float)
        #otu_table = convert_otu_table_relative(otu_table)
        otu_table = self.otu_table.normObservationBySample()
        summary, header = make_summary(otu_table, 3, upper_percentage, lower_percentage)
        self.assertEqual(header, ['Taxon','s1','s2','s3','s4'])
        self.assertEqual(summary[0][0], ('Root','Bacteria','Actinobacteria'))
        self.assertFloatEqual(summary[0][1:], [1.0 / 3, 0.0, 0.5, 0.8])
        self.assertEqual(summary[1][0], ('Root','Bacteria','Firmicutes'))
        self.assertFloatEqual(summary[1][1:], [1.0 / 3, 0.6, 0.25, 0.2])
        self.assertEqual(summary[2][0], ('Root','Bacteria','Other'))
        self.assertFloatEqual(summary[2][1:], [1.0 / 3, 0.4, 0.25, 0.0])
        
        ##
        # testing lower triming 
        lower_percentage, upper_percentage = 0.3, None
        summary, header = make_summary(otu_table, 3, upper_percentage, lower_percentage)
        self.assertEqual(summary[0][0], ('Root','Bacteria','Other'))
        self.assertFloatEqual(summary[0][1:], [1.0 / 3, 0.4, 0.25, 0.0])
        
        ##
        # testing upper triming 
        lower_percentage, upper_percentage = None, 0.4
        summary, header = make_summary(otu_table, 3, upper_percentage, lower_percentage)
        self.assertEqual(summary[0][0], ('Root','Bacteria','Actinobacteria'))
        self.assertFloatEqual(summary[0][1:], [1.0 / 3, 0.0, 0.5, 0.8])
        

    def test_add_summary_category_mapping(self):
        """make_new_summary_file works
        """
        #otu_table = parse_otu_table(self.otu_table, int)
        #otu_table = parse_biom_table(self.otu_table)
        mapping, header, comments = parse_mapping_file(self.mapping)
        summary, taxon_order = add_summary_mapping(self.otu_table, mapping, 3)
        self.assertEqual(taxon_order, [('Root','Bacteria','Actinobacteria'),
                                       ('Root','Bacteria','Firmicutes'),
                                       ('Root','Bacteria','Other')])
        self.assertEqual(summary, {'s1':[1,1,1],
                                   's2':[0,3,2],
                                   's3':[2,1,1],
                                   's4':[4,1,0]})

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
