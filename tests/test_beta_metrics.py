#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.0-rc2"
__maintainer__ = "justin kuczynski"
__email__ = "justinak@gmail.com"


"""Contains tests for beta_metrics  functions."""
import os.path
import numpy
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from qiime.parse import make_envs_dict
from qiime.beta_metrics import (
    _reorder_unifrac_res,
    make_unifrac_metric,
    make_unifrac_row_metric)
from qiime.parse import parse_newick
from cogent.core.tree import PhyloNode
from cogent.maths.unifrac.fast_tree import (unifrac)
import warnings


class FunctionTests(TestCase):

    def setUp(self):
        self.l19_data = numpy.array([
            [7, 1, 0, 0, 0, 0, 0, 0, 0],
            [4, 2, 0, 0, 0, 1, 0, 0, 0],
            [2, 4, 0, 0, 0, 1, 0, 0, 0],
            [1, 7, 0, 0, 0, 0, 0, 0, 0],
            [0, 8, 0, 0, 0, 0, 0, 0, 0],
            [0, 7, 1, 0, 0, 0, 0, 0, 0],
            [0, 4, 2, 0, 0, 0, 2, 0, 0],
            [0, 2, 4, 0, 0, 0, 1, 0, 0],
            [0, 1, 7, 0, 0, 0, 0, 0, 0],
            [0, 0, 8, 0, 0, 0, 0, 0, 0],
            [0, 0, 7, 1, 0, 0, 0, 0, 0],
            [0, 0, 4, 2, 0, 0, 0, 3, 0],
            [0, 0, 2, 4, 0, 0, 0, 1, 0],
            [0, 0, 1, 7, 0, 0, 0, 0, 0],
            [0, 0, 0, 8, 0, 0, 0, 0, 0],
            [0, 0, 0, 7, 1, 0, 0, 0, 0],
            [0, 0, 0, 4, 2, 0, 0, 0, 4],
            [0, 0, 0, 2, 4, 0, 0, 0, 1],
            [0, 0, 0, 1, 7, 0, 0, 0, 0]
        ])
        self.l19_sample_names = [
            'sam1', 'sam2', 'sam3', 'sam4', 'sam5', 'sam6',
            'sam7', 'sam8', 'sam9', 'sam_middle', 'sam11', 'sam12', 'sam13',
            'sam14', 'sam15', 'sam16', 'sam17', 'sam18', 'sam19']
        self.l19_taxon_names = ['tax1', 'tax2', 'tax3', 'tax4', 'endbigtaxon',
                                'tax6', 'tax7', 'tax8', 'tax9']
        self.l19_treestr = '((((tax7:0.1,tax3:0.2):.98,tax8:.3, tax4:.3):.4, ' +\
            '((tax1:0.3, tax6:.09):0.43,tax2:0.4):0.5):.2,' +\
            '(tax9:0.3, endbigtaxon:.08));'

    def test_reorder_unifrac_res(self):
        """ reorder_unifrac_res should correctly reorder a misordered 3x3 matrix
        """
        mtx = numpy.array([[1, 2, 3],
                           [4, 5, 6],
                           [7, 8, 9]], 'float')
        unifrac_mtx = numpy.array([[1, 3, 2],
                                   [7, 9, 8],
                                   [4, 6, 5]], 'float')
        sample_names = ['yo', "it's", "samples"]
        unifrac_sample_names = ['yo', "samples", "it's"]
        reordered_mtx = _reorder_unifrac_res(
            [unifrac_mtx, unifrac_sample_names],
            sample_names)
        assert_almost_equal(reordered_mtx, mtx)

    def test_make_unifrac_metric(self):
        """ exercise of the unweighted unifrac metric should not throw errors"""
        tree = parse_newick(self.l19_treestr, PhyloNode)
        unif = make_unifrac_metric(False, unifrac, True)
        res = unif(self.l19_data, self.l19_taxon_names, tree,
                   self.l19_sample_names)
        envs = make_envs_dict(self.l19_data, self.l19_sample_names,
                              self.l19_taxon_names)
        unifrac_mat, unifrac_names = fast_unifrac(tree, envs,
                                                  modes=['distance_matrix'])['distance_matrix']
        assert_almost_equal(res, _reorder_unifrac_res([unifrac_mat,
                                                         unifrac_names], self.l19_sample_names))
        self.assertEqual(res[0, 0], 0)
        self.assertEqual(res[0, 3], 0.0)
        self.assertNotEqual(res[0, 1], 1.0)

    def test_make_unifrac_metric2(self):
        """ samples with no seqs, and identical samples, should behave correctly
        """
        tree = parse_newick(self.l19_treestr, PhyloNode)
        unif = make_unifrac_metric(False, unifrac, True)
        otu_data = numpy.array([
            [0, 0, 0, 0, 0, 0, 0, 0, 0],  # sam1 zeros
            [4, 2, 0, 0, 0, 1, 0, 0, 0],
            [2, 4, 0, 0, 0, 1, 0, 0, 0],
            [1, 7, 0, 0, 0, 0, 0, 0, 0],
            [0, 8, 0, 0, 0, 0, 0, 0, 0],
            [0, 7, 1, 0, 0, 0, 0, 0, 0],
            [0, 4, 2, 0, 0, 0, 2, 0, 0],
            [0, 2, 4, 0, 0, 0, 1, 0, 0],
            [0, 1, 7, 0, 0, 0, 0, 0, 0],
            [0, 0, 8, 0, 0, 0, 0, 0, 0],
            [0, 0, 7, 1, 0, 0, 0, 0, 0],
            [0, 0, 4, 2, 0, 0, 0, 3, 0],
            [0, 0, 2, 4, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],  # sam14 zeros
            [0, 0, 0, 8, 0, 0, 0, 0, 0],
            [0, 0, 2, 4, 0, 0, 0, 1, 0],  # sam 16 now like sam 13
            [0, 0, 0, 4, 2, 0, 0, 0, 4],
            [0, 0, 0, 2, 4, 0, 0, 0, 1],
            [0, 0, 0, 1, 7, 0, 0, 0, 0]
        ])
        warnings.filterwarnings('ignore')
        res = unif(otu_data, self.l19_taxon_names, tree,
                   self.l19_sample_names)
        envs = make_envs_dict(self.l19_data, self.l19_sample_names,
                              self.l19_taxon_names)
        self.assertEqual(res[0, 0], 0)
        self.assertEqual(res[0, 13], 0.0)
        self.assertEqual(res[12, 15], 0.0)
        self.assertEqual(res[0, 1], 1.0)
        warnings.resetwarnings()

    def test_make_unifrac_metric3(self):
        treestr = '((((tax7:0.1):.98,tax8:.3, tax4:.3):.4, ' +\
            '((tax6:.09):0.43):0.5):.2,' +\
            '(tax9:0.3, endbigtaxon:.08));'  # taxa 1,2,3 removed
        tree = parse_newick(treestr, PhyloNode)

        otu_data = numpy.array([
            [7, 1, 0, 0, 0, 0, 0, 0, 0],  # 1 now zeros
            [4, 2, 0, 0, 0, 1, 0, 0, 0],
            [2, 4, 0, 0, 0, 1, 0, 0, 0],
            [1, 7, 0, 0, 0, 0, 0, 0, 0],  # 4 now zeros
            [0, 8, 0, 0, 0, 0, 0, 0, 0],
            [0, 7, 1, 0, 0, 0, 0, 0, 0],
            [0, 4, 2, 0, 0, 0, 2, 0, 0],
            [0, 2, 4, 0, 0, 0, 1, 0, 0],
            [0, 1, 7, 0, 0, 0, 0, 0, 0],
            [0, 0, 8, 0, 0, 0, 0, 0, 0],
            [0, 0, 7, 1, 0, 0, 0, 0, 0],
            [0, 0, 4, 2, 0, 0, 0, 3, 0],
            [0, 0, 2, 4, 0, 0, 0, 1, 0],
            [0, 0, 1, 7, 0, 0, 0, 0, 0],
            [0, 0, 0, 8, 0, 0, 0, 0, 0],
            [0, 0, 0, 7, 1, 0, 0, 0, 0],
            [0, 0, 0, 4, 2, 0, 0, 0, 4],
            [0, 0, 0, 2, 4, 0, 0, 0, 1],
            [0, 0, 0, 1, 7, 0, 0, 0, 0]
        ])

        unif = make_unifrac_metric(False, unifrac, True)
        warnings.filterwarnings('ignore')
        res = unif(otu_data, self.l19_taxon_names, tree,
                   self.l19_sample_names)
        warnings.resetwarnings()
        envs = make_envs_dict(self.l19_data, self.l19_sample_names,
                              self.l19_taxon_names)
        self.assertEqual(res[0, 0], 0)
        self.assertEqual(res[0, 3], 0.0)
        self.assertEqual(res[0, 1], 1.0)

    def test_make_unifrac_row_metric3(self):
        treestr = '((((tax7:0.1):.98,tax8:.3, tax4:.3):.4, ' +\
            '((tax6:.09):0.43):0.5):.2,' +\
            '(tax9:0.3, endbigtaxon:.08));'  # taxa 1,2,3 removed
        tree = parse_newick(treestr, PhyloNode)

        otu_data = numpy.array([
            [7, 1, 0, 0, 0, 0, 0, 0, 0],  # 1 now zeros
            [4, 2, 0, 0, 0, 1, 0, 0, 0],
            [2, 4, 0, 0, 0, 1, 0, 0, 0],
            [1, 7, 0, 0, 0, 0, 0, 0, 0],  # 4 now zeros
            [0, 8, 0, 0, 0, 0, 0, 0, 0],
            [0, 7, 1, 0, 0, 0, 0, 0, 0],
            [0, 4, 2, 0, 0, 0, 2, 0, 0],
            [0, 2, 4, 0, 0, 0, 1, 0, 0],
            [0, 1, 7, 0, 0, 0, 0, 0, 0],
            [0, 0, 8, 0, 0, 0, 0, 0, 0],
            [0, 0, 7, 1, 0, 0, 0, 0, 0],
            [0, 0, 4, 2, 0, 0, 0, 3, 0],
            [0, 0, 2, 4, 0, 0, 0, 1, 0],
            [0, 0, 1, 7, 0, 0, 0, 0, 0],
            [0, 0, 0, 8, 0, 0, 0, 0, 0],
            [0, 0, 0, 7, 1, 0, 0, 0, 0],
            [0, 0, 0, 4, 2, 0, 0, 0, 4],
            [0, 0, 0, 2, 4, 0, 0, 0, 1],
            [0, 0, 0, 1, 7, 0, 0, 0, 0]
        ])

        unif = make_unifrac_metric(False, unifrac, True)
        warnings.filterwarnings('ignore')
        res = unif(otu_data, self.l19_taxon_names, tree,
                   self.l19_sample_names)
        warnings.resetwarnings()
        envs = make_envs_dict(self.l19_data, self.l19_sample_names,
                              self.l19_taxon_names)
        self.assertEqual(res[0, 0], 0)
        self.assertEqual(res[0, 3], 0.0)
        self.assertEqual(res[0, 1], 1.0)

        warnings.filterwarnings('ignore')
        unif_row = make_unifrac_row_metric(False, unifrac, True)
        for i, sam_name in enumerate(self.l19_sample_names):
            if i in [0, 3, 4, 5, 8, 9]:
                continue
            # these have no data and are warned "meaningless".
            # I Would prefer if they matched res anyway though
            res_row = unif_row(otu_data, self.l19_taxon_names, tree,
                               self.l19_sample_names, sam_name)
            for j in range(len(self.l19_sample_names)):
                if j in [0, 3, 4, 5, 8, 9]:
                    continue  # ok if meaningless number in zero sample
                self.assertAlmostEqual(res_row[j], res[i, j])
        warnings.resetwarnings()

    def test_make_unifrac_row_metric2(self):
        """ samples with no seqs, and identical samples, should behave correctly
        """
        tree = parse_newick(self.l19_treestr, PhyloNode)
        unif = make_unifrac_metric(False, unifrac, True)
        otu_data = numpy.array([
            [0, 0, 0, 0, 0, 0, 0, 0, 0],  # sam1 zeros
            [4, 2, 0, 0, 0, 1, 0, 0, 0],
            [2, 4, 0, 0, 0, 1, 0, 0, 0],
            [1, 7, 0, 0, 0, 0, 0, 0, 0],
            [0, 8, 0, 0, 0, 0, 0, 0, 0],
            [0, 7, 1, 0, 0, 0, 0, 0, 0],
            [0, 4, 2, 0, 0, 0, 2, 0, 0],
            [0, 2, 4, 0, 0, 0, 1, 0, 0],
            [0, 1, 7, 0, 0, 0, 0, 0, 0],
            [0, 0, 8, 0, 0, 0, 0, 0, 0],
            [0, 0, 7, 1, 0, 0, 0, 0, 0],
            [0, 0, 4, 2, 0, 0, 0, 3, 0],
            [0, 0, 2, 4, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0],  # sam14 zeros
            [0, 0, 0, 8, 0, 0, 0, 0, 0],
            [0, 0, 2, 4, 0, 0, 0, 1, 0],  # sam 16 now like sam 13
            [0, 0, 0, 4, 2, 0, 0, 0, 4],
            [0, 0, 0, 2, 4, 0, 0, 0, 1],
            [0, 0, 0, 1, 7, 0, 0, 0, 0]
        ])
        warnings.filterwarnings('ignore')
        res = unif(otu_data, self.l19_taxon_names, tree,
                   self.l19_sample_names)
        envs = make_envs_dict(self.l19_data, self.l19_sample_names,
                              self.l19_taxon_names)
        self.assertEqual(res[0, 0], 0)
        self.assertEqual(res[0, 13], 0.0)
        self.assertEqual(res[12, 15], 0.0)
        self.assertEqual(res[0, 1], 1.0)
        warnings.resetwarnings()

        warnings.filterwarnings('ignore')
        unif_row = make_unifrac_row_metric(False, unifrac, True)
        for i, sam_name in enumerate(self.l19_sample_names):
            if i in [0]:
                continue
            # these have no data and are warned "meaningless".
            # I Would prefer if they matched res anyway though
            res_row = unif_row(otu_data, self.l19_taxon_names, tree,
                               self.l19_sample_names, sam_name)
            for j in range(len((self.l19_sample_names))):
                if j in [0]:
                    continue  # ok if meaningless number in zero sample
                self.assertEqual(res_row[j], res[i, j])
        warnings.resetwarnings()

# run tests if called from command line
if __name__ == '__main__':
    main()
