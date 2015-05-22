#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""Contains tests for performing alpha diversity analyses within each sample."""

from os import makedirs, close
from tempfile import mkstemp
from unittest import TestCase, main

from biom.table import Table
from cogent.maths.unifrac.fast_unifrac import PD_whole_tree
from numpy import array
from numpy.testing import assert_almost_equal
from skbio.diversity.alpha import observed_otus, osd
from skbio.util import remove_files

from qiime.alpha_diversity import (AlphaDiversityCalc, AlphaDiversityCalcs,
                                   single_file_cup)
from qiime.parse import parse_newick
from qiime.util import get_qiime_temp_dir, write_biom_table


class AlphaDiversitySharedSetUpTests(TestCase):

    def setUp(self):
        """Define some test data."""
        self.tmp_dir = get_qiime_temp_dir()

        self.otu_table1 = Table(data=array([[2, 0, 0, 1],
                                            [1, 1, 1, 1],
                                            [0, 0, 0, 0]]).T,
                                           sample_ids=list('XYZ'),
                                           observation_ids=list('abcd'))
        fd, self.otu_table1_fp = mkstemp(dir=self.tmp_dir,
                                              prefix='alpha_diversity_tests',
                                              suffix='.biom')
        close(fd)
        write_biom_table(self.otu_table1, self.otu_table1_fp)

        self.otu_table2 = Table(data=array([[2, 0, 0, 1],
                                                   [1, 1, 1, 1],
                                                   [0, 0, 0, 0]]).T,
                                        sample_ids=list('XYZ'),
                                        observation_ids=['a', 'b', 'c', 'd_'])
        fd, self.otu_table2_fp = mkstemp(dir=self.tmp_dir,
                                              prefix='alpha_diversity_tests',
                                              suffix='.biom')
        close(fd)
        write_biom_table(self.otu_table2, self.otu_table2_fp)

        self.single_sample_otu_table = Table(
            data=array([[2, 0, 0, 1]]).T,
            sample_ids=list('X'),
            observation_ids=list(
                'abcd'))
        fd, self.single_sample_otu_table_fp = mkstemp(
            dir=self.tmp_dir,
            prefix='alpha_diversity_tests',
            suffix='.biom')
        close(fd)
        write_biom_table(self.single_sample_otu_table,
                         self.single_sample_otu_table_fp)

        self.tree1 = parse_newick('((a:2,b:3):2,(c:1,d:2):7);')
        self.tree2 = parse_newick("((a:2,'b':3):2,(c:1,'d_':2):7);")

        self.files_to_remove = [self.otu_table1_fp, self.otu_table2_fp,
                                self.single_sample_otu_table_fp]

    def tearDown(self):
        """ """
        remove_files(self.files_to_remove)


class AlphaDiversityCalcTests(AlphaDiversitySharedSetUpTests):

    """Tests of the AlphaDiversityCalc class"""

    def test_init(self):
        """AlphaDiversity __init__ should store metric, name, params"""
        c = AlphaDiversityCalc(observed_otus)
        self.assertEqual(c.Metric, observed_otus)
        self.assertEqual(c.Params, {})

    def test_call(self):
        """AlphaDiversityCalc __call__ should call metric on data
        and return correct result"""
        c = AlphaDiversityCalc(observed_otus)
        assert_almost_equal(c(data_path=self.otu_table1_fp), [2, 4, 0])

    def test_multi_return(self):
        """AlphaDiversityCalc __call__ should call metric on data
        and return correct result for metric fn that returns a len 3 tuple
        """
        c = AlphaDiversityCalc(osd)
        res = c(data_path=self.otu_table1_fp)
        assert_almost_equal(res, array([[2, 1, 1],
                            [4, 4, 0],
                            [0, 0, 0]]))

    def test_1sample(self):
        """ should work if only testing one sample as well"""
        c = AlphaDiversityCalc(observed_otus)
        self.assertEqual(c(data_path=self.single_sample_otu_table_fp), [2])

    def test_call_phylogenetic(self):
        """AlphaDiversityCalc __call__ should call metric on phylo data
        and return correct values"""
        c = AlphaDiversityCalc(metric=PD_whole_tree,
                               is_phylogenetic=True)
        assert_almost_equal(
            c(data_path=self.otu_table1_fp, tree_path=self.tree1,
              taxon_names=self.otu_table1.ids(axis='observation'),
              sample_names=self.otu_table1.ids()), [13, 17, 0])

    def test_call_phylogenetic_escaped_names(self):
        """AlphaDiversityCalc __call__ should call metric on phylo data
        and return correct values"""

        c = AlphaDiversityCalc(metric=PD_whole_tree, is_phylogenetic=True)
        expected = [13., 17., 0.]

        non_escaped_result = c(
            data_path=self.otu_table1_fp, tree_path=self.tree1,
            taxon_names=self.otu_table1.ids(axis='observation'),
            sample_names=self.otu_table1.ids())

        escaped_result = c(data_path=self.otu_table2_fp,
                           tree_path=self.tree2,
                           taxon_names=self.otu_table2.ids(axis='observation'),
                           sample_names=self.otu_table2.ids())

        assert_almost_equal(non_escaped_result, expected)
        assert_almost_equal(escaped_result, expected)
        assert_almost_equal(non_escaped_result, escaped_result)


class AlphaDiversityCalcsTests(AlphaDiversitySharedSetUpTests):

    """Tests of the AlphaDiversityCalcs class"""

    def test1(self):
        """ checks that output from AlphaDiversityCalcs is the right shape
        when run on phylo, multiple return value nonphylo, and another nonphylo
        """
        calc1 = AlphaDiversityCalc(metric=observed_otus)
        calc2 = AlphaDiversityCalc(metric=PD_whole_tree,
                                   is_phylogenetic=True)
        calc3 = AlphaDiversityCalc(metric=osd)
        calcs = AlphaDiversityCalcs([calc1, calc2, calc3])
        results = calcs(data_path=self.otu_table1_fp,
                        tree_path=self.tree1,
                        result_path=None, log_path=None)
        self.assertEqual(results[0].shape, (3, 5))
        self.assertEqual(len(results[1]), 3)
        self.assertEqual(len(results[2]), 5)


class SingleFileCUPTests(TestCase):
    def setUp(self):
        self.files_to_remove = []

        fd, self.tmp_file = mkstemp(suffix="test_single_file_cup.biom")
        close(fd)
        self.files_to_remove.append(self.tmp_file)

        fd, self.tmp_outfile = mkstemp(suffix="test_single_file_cup.txt")
        close(fd)
        self.files_to_remove.append(self.tmp_outfile)

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_single_file_cup_string(self):
        """Returns matrix with estimates using metrics string."""
        # convert_biom using otu_table w/o leading #
        bt_string = (
            '{"rows": [{"id": "1", "metadata": null}, {"id": "2",'
            '"metadata": null}, {"id": "3", "metadata": null}, {"id": "4", '
            '"metadata": null}, {"id": "5", "metadata": null}], "format": '
            '"Biological Observation Matrix 0.9.1-dev", "data": [[0, 0, 3.0], '
            '[0, 1, 4.0], [1, 0, 2.0], [1, 1, 5.0], [2, 0, 1.0], [2, 1, 2.0], '
            '[3, 1, 4.0], [4, 0, 1.0]], "columns": [{"id": "S1", "metadata": '
            'null}, {"id": "S2", "metadata": null}], "generated_by": '
            '"BIOM-Format 0.9.1-dev", "matrix_type": "sparse", "shape": '
            '[5, 2], "format_url": "http://biom-format.org", "date": '
            '"2012-05-04T09:28:28.247809", "type": "OTU table", "id": null, '
            '"matrix_element_type": "float"}')

        with open(self.tmp_file, 'w') as fh:
            fh.write(bt_string)

        single_file_cup(self.tmp_file, 'lladser_pe,lladser_ci',
                        self.tmp_outfile, r=4)

        # Not much testing here, just make sure we get back a (formatted)
        # matrix with the right dimensions
        with open(self.tmp_outfile, 'U') as out_f:
            observed = out_f.readlines()
        self.assertEqual(len(observed), 3)
        self.assertEqual(len(observed[1].split('\t')), 4)

    def test_single_file_cup_list(self):
        """Returns matrix with estimates using metrics list."""
        # convert_biom using otu_table w/o leading #
        bt_string = (
            '{"rows": [{"id": "1", "metadata": null}], "format": "Biological '
            'Observation Matrix 0.9.1-dev", "data": [[0, 0, 3.0]], "columns": '
            '[{"id": "S1", "metadata": null}], "generated_by": '
            '"BIOM-Format 0.9.1-dev", "matrix_type": "sparse", "shape": '
            '[1, 1], "format_url": "http://biom-format.org", "date": '
            '"2012-05-04T09:36:57.500673", "type": "OTU table", "id": null, '
            '"matrix_element_type": "float"}')

        with open(self.tmp_file, 'w') as fh:
            fh.write(bt_string)

        single_file_cup(self.tmp_file, ['lladser_pe', 'lladser_ci'],
                        self.tmp_outfile, r=4)

        with open(self.tmp_outfile, 'U') as out_f:
            observed = out_f.readlines()

        expected = ["\tlladser_pe\tlladser_lower_bound\tlladser_upper_bound\n",
                    "S1\tnan\tnan\tnan\n"]
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    main()
