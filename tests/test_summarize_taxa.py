#!/usr/bin/env python
from __future__ import division

"""Tests of code for summarizing taxa in an OTU table"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Daniel McDonald", "Antonio Gonzalez Pena",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_mapping_file
from qiime.util import convert_otu_table_relative
from numpy import array
from biom.exception import TableException
from biom.table import SparseOTUTable, SparseTaxonTable, table_factory
from biom.parse import parse_biom_table
from qiime.summarize_taxa import (_make_collapse_fn, make_summary,
                                  add_summary_mapping)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        # #OTU ID s1 s2 s3 s4 Consensus Lineage
        # 0 1 0 2 4 Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
        # 1 1 2 0 1 Root;Bacteria;Firmicutes;"Clostridia"
        # 2 0 1 1 0 Root;Bacteria;Firmicutes;"Clostridia"
        # 3 1 2 1 0 Root;Bacteria

        otu_table_vals = array([[1,0,2,4],
                                [1,2,0,1],
                                [0,1,1,0],
                                [1,2,1,0]])
        sample_ids = ['s1', 's2', 's3', 's4']
        obs_ids = ['0', '1', '2', '3']

        md_as_list = [{"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]},
                      {"taxonomy": ["Root", "Bacteria", "Firmicutes", "\"Clostridia\""]},
                      {"taxonomy": ["Root", "Bacteria", "Firmicutes", "\"Clostridia\""]},
                      {"taxonomy": ["Root", "Bacteria"]}]

        md_as_string = [{"taxonomy": "Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae"},
                        {"taxonomy": "Root;Bacteria;Firmicutes;\"Clostridia\""},
                        {"taxonomy": "Root;Bacteria;Firmicutes;\"Clostridia\""},
                        {"taxonomy": "Root;Bacteria"}]

        # Mixed 1-1 and 1-M metadata, in various supported formats.
        one_to_many_md = [{"taxonomy": [['a', 'b', 'c'],
                                        ['a', 'b']]},
                          {"taxonomy": ['a', 'b', 'c']},
                          {"taxonomy": [['a', 'bb', 'c', 'd']]},
                          {"taxonomy": [['a', 'bb'], ['b']]}]

        self.otu_table = table_factory(otu_table_vals, sample_ids, obs_ids,
                                       None, md_as_list)

        self.otu_table_rel = self.otu_table.normObservationBySample()

        self.otu_table_md_as_string = table_factory(otu_table_vals, sample_ids,
                                                    obs_ids, None,
                                                    md_as_string)

        self.minimal_table = table_factory(otu_table_vals, sample_ids, obs_ids,
                                           None, None)

        self.otu_table_one_to_many = table_factory(otu_table_vals, sample_ids,
                                                   obs_ids, None,
                                                   one_to_many_md)

        self.mapping="""#SampleID\tBarcodeSequence\tTreatment\tDescription
#Test mapping file
s1\tAAAA\tControl\tControl mouse, I.D. 354
s2\tGGGG\tControl\tControl mouse, I.D. 355
s3\tCCCC\tExp\tDisease mouse, I.D. 356
s4\tTTTT\tExp\tDisease mouse, I.D. 357""".split('\n')

    def test_make_summary(self):
        """make_summary works"""
        # level 2
        exp_data = array([3.0, 5.0, 4.0, 5.0])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root;Bacteria'])
        obs = make_summary(self.otu_table, 2, absolute_abundance=True)
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), SparseTaxonTable)

        # level 3
        exp_data = array([[1.0, 0.0, 2.0, 4.0],
                          [1.0, 3.0, 1.0, 1.0],
                          [1.0, 2.0, 1.0, 0.0]])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root;Bacteria;Actinobacteria',
                             'Root;Bacteria;Firmicutes',
                             'Root;Bacteria;Other'])
        obs = make_summary(self.otu_table, 3, absolute_abundance=True)
        self.assertEqual(obs, exp)

        # level 4
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root;Bacteria;Actinobacteria;Actinobacteria',
                             'Root;Bacteria;Firmicutes;"Clostridia"',
                             'Root;Bacteria;Other;Other'])
        obs = make_summary(self.otu_table, 4, absolute_abundance=True)
        self.assertEqual(obs, exp)

        # md_as_string=True
        obs = make_summary(self.otu_table_md_as_string, 4,
                           absolute_abundance=True, md_as_string=True)
        self.assertEqual(obs, exp)

        # custom delimiter
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root>Bacteria>Actinobacteria>Actinobacteria',
                             'Root>Bacteria>Firmicutes>"Clostridia"',
                             'Root>Bacteria>Other>Other'])
        obs = make_summary(self.otu_table, 4, absolute_abundance=True,
                           delimiter='>')
        self.assertEqual(obs, exp)

        # custom constructor
        obs = make_summary(self.otu_table, 4, absolute_abundance=True,
                           delimiter='>', constructor=SparseOTUTable)
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), SparseOTUTable)

        # absolute_abudance=False
        exp_data = array([1.0, 1.0, 1.0, 1.0])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root;Bacteria'])
        obs = make_summary(self.otu_table, 2, absolute_abundance=False)
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), SparseTaxonTable)

    def test_make_summary_invalid_input(self):
        """make_summary handles invalid input"""
        # No metadata.
        with self.assertRaises(ValueError):
            make_summary(self.minimal_table, 2)

        # Wrong metadata key.
        with self.assertRaises(KeyError):
            make_summary(self.otu_table, 2, md_identifier='foo')

    def test_make_summary_relative_abundances(self):
        """make_summary works with relative abundances"""
        exp_data = array([[1.0 / 3, 0.0, 0.5, 0.8],
                          [1.0 / 3, 0.6, 0.25, 0.2],
                          [1.0 / 3, 0.4, 0.25, 0.0]])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root;Bacteria;Actinobacteria',
                             'Root;Bacteria;Firmicutes',
                             'Root;Bacteria;Other'])

        obs = make_summary(self.otu_table_rel, 3)

        # Can't use __eq__ here because of floating point error.
        self.assertEqual(obs.SampleIds, exp.SampleIds)
        self.assertEqual(obs.ObservationIds, exp.ObservationIds)
        self.assertFloatEqual(obs.sampleData('s1'), exp.sampleData('s1'))
        self.assertFloatEqual(obs.sampleData('s2'), exp.sampleData('s2'))
        self.assertFloatEqual(obs.sampleData('s3'), exp.sampleData('s3'))
        self.assertFloatEqual(obs.sampleData('s4'), exp.sampleData('s4'))

    def test_make_summary_trimming(self):
        """make_summary correctly trims taxa based on abundance"""
        # testing lower trimming 
        exp_data = array([[1.0 / 3, 0.4, 0.25, 0.0]])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root;Bacteria;Other'])

        obs = make_summary(self.otu_table_rel, 3, absolute_abundance=True,
                           lower_percentage=0.3)

        self.assertEqual(obs.SampleIds, exp.SampleIds)
        self.assertEqual(obs.ObservationIds, exp.ObservationIds)
        self.assertFloatEqual(obs.sampleData('s1'), exp.sampleData('s1'))
        self.assertFloatEqual(obs.sampleData('s2'), exp.sampleData('s2'))
        self.assertFloatEqual(obs.sampleData('s3'), exp.sampleData('s3'))
        self.assertFloatEqual(obs.sampleData('s4'), exp.sampleData('s4'))

        # testing upper trimming 
        exp_data = array([[1.0 / 3, 0.0, 0.5, 0.8]])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root;Bacteria;Actinobacteria'])

        obs = make_summary(self.otu_table_rel, 3, absolute_abundance=True,
                           upper_percentage=0.4)

        self.assertEqual(obs.SampleIds, exp.SampleIds)
        self.assertEqual(obs.ObservationIds, exp.ObservationIds)
        self.assertFloatEqual(obs.sampleData('s1'), exp.sampleData('s1'))
        self.assertFloatEqual(obs.sampleData('s2'), exp.sampleData('s2'))
        self.assertFloatEqual(obs.sampleData('s3'), exp.sampleData('s3'))
        self.assertFloatEqual(obs.sampleData('s4'), exp.sampleData('s4'))

        # test lower and upper trimming 
        exp_data = array([[1.0 / 3, 0.6, 0.25, 0.2]])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['Root;Bacteria;Firmicutes'])

        obs = make_summary(self.otu_table_rel, 3, absolute_abundance=True,
                           upper_percentage=0.3, lower_percentage=0.4)

        self.assertEqual(obs.SampleIds, exp.SampleIds)
        self.assertEqual(obs.ObservationIds, exp.ObservationIds)
        self.assertFloatEqual(obs.sampleData('s1'), exp.sampleData('s1'))
        self.assertFloatEqual(obs.sampleData('s2'), exp.sampleData('s2'))
        self.assertFloatEqual(obs.sampleData('s3'), exp.sampleData('s3'))
        self.assertFloatEqual(obs.sampleData('s4'), exp.sampleData('s4'))

        # test trimming everything out
        with self.assertRaises(TableException):
            make_summary(self.otu_table_rel, 3, upper_percentage=0.2,
                         lower_percentage=0.2)

    def test_make_summary_one_to_many(self):
        """make_summary works with one-to-many obs-md relationship"""
        # one_to_many='first'
        exp_data = array([[2.0, 2.0, 2.0, 5.0],
                          [1.0, 2.0, 1.0, 0.0],
                          [0.0, 1.0, 1.0, 0.0]])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['a;b;c', 'a;bb;Other', 'a;bb;c'])

        obs = make_summary(self.otu_table_one_to_many, 3,
                           absolute_abundance=True, one_to_many='first')
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), SparseTaxonTable)

        # one_to_many='add'
        exp_data = array([[1.0, 0.0, 2.0, 4.0],
                          [2.0, 2.0, 2.0, 5.0],
                          [1.0, 2.0, 1.0, 0.0],
                          [0.0, 1.0, 1.0, 0.0],
                          [1.0, 2.0, 1.0, 0.0]])
        exp = table_factory(exp_data, ['s1', 's2', 's3', 's4'],
                            ['a;b;Other', 'a;b;c', 'a;bb;Other', 'a;bb;c',
                             'b;Other;Other'])

        obs = make_summary(self.otu_table_one_to_many, 3,
                           absolute_abundance=True, one_to_many='add')
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), SparseTaxonTable)

    def test_add_summary_mapping(self):
        """add_summary_mapping works"""
        mapping, header, comments = parse_mapping_file(self.mapping)
        summary, taxon_order = add_summary_mapping(self.otu_table, mapping, 3,
                                                   absolute_abundance=True,
                                                   delimiter='FOO')
        self.assertEqual(taxon_order, ('RootFOOBacteriaFOOActinobacteria',
                                       'RootFOOBacteriaFOOFirmicutes',
                                       'RootFOOBacteriaFOOOther'))
        self.assertEqual(summary, {'s1':[1,1,1],
                                   's2':[0,3,2],
                                   's3':[2,1,1],
                                   's4':[4,1,0]})

    def test_make_collapse_fn_invalid_input(self):
        """_make_collapse_fn correctly handles invalid input"""
        with self.assertRaises(ValueError):
            _make_collapse_fn(1, one_to_many='foo')


#run unit tests if run from command-line
if __name__ == '__main__':
    main()
