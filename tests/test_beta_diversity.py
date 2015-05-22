#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Jose Antonio Navas Molina",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


"""Contains tests for performing beta diversity analyses."""

import warnings
import os
import shutil
from tempfile import mkstemp, mkdtemp
from unittest import TestCase, main

from biom.table import Table
import numpy as np
import numpy.testing as npt
from skbio.util import remove_files
from cogent.core.tree import PhyloNode
from cogent.maths.distance_transform import dist_chisq

from qiime.util import get_qiime_temp_dir, write_biom_table
from qiime.parse import parse_newick, parse_distmat, parse_matrix
from qiime.beta_diversity import BetaDiversityCalc, single_file_beta,\
    list_known_nonphylogenetic_metrics, list_known_phylogenetic_metrics
from qiime.beta_metrics import dist_unweighted_unifrac


class BetaDiversityCalcTests(TestCase):

    """Tests of the BetaDiversityCalc class"""

    def setUp(self):
        self.tmp_dir = get_qiime_temp_dir()

        self.l19_data = np.array([
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
        self.l19_taxon_names_w_underscore = ['ta_x1', 'tax2', 'tax3', 'tax4',
                                             'endbigtaxon', 'tax6', 'tax7',
                                             'tax8', 'tax9']

        l19 = Table(self.l19_data.T, self.l19_taxon_names,
                    self.l19_sample_names)
        fd, self.l19_fp = mkstemp(dir=self.tmp_dir,
                                prefix='test_bdiv_otu_table', suffix='.blom')
        os.close(fd)
        write_biom_table(l19, self.l19_fp)

        l19_w_underscore = Table(self.l19_data.T,
                                 self.l19_taxon_names_w_underscore,
                                 self.l19_sample_names)
        fd, self.l19_w_underscore_fp = mkstemp(dir=self.tmp_dir,
                                               prefix='test_bdiv_otu_table',
                                               suffix='.blom')
        os.close(fd)
        write_biom_table(l19_w_underscore, self.l19_w_underscore_fp)

        self.l19_tree_str = '((((tax7:0.1,tax3:0.2):.98,tax8:.3, tax4:.3):.4,\
 ((tax1:0.3, tax6:.09):0.43,tax2:0.4):0.5):.2, (tax9:0.3, endbigtaxon:.08));'
        self.l19_tree = parse_newick(self.l19_tree_str, PhyloNode)

        self.files_to_remove = [self.l19_fp, self.l19_w_underscore_fp]
        self.folders_to_remove = []

    def tearDown(self):
        remove_files(self.files_to_remove)
        for folder in self.folders_to_remove:
            shutil.rmtree(folder)

    def test_l19_chi(self):
        """beta calc run should return same values as directly calling metric"""
        beta_calc_chisq = BetaDiversityCalc(dist_chisq, 'chi square', False)
        matrix, labels = beta_calc_chisq(data_path=self.l19_fp, tree_path=None)
        self.assertEqual(labels, self.l19_sample_names)
        npt.assert_almost_equal(matrix, dist_chisq(self.l19_data))

    def test_l19_unifrac(self):
        """beta calc run should also work for phylo metric"""
        beta_calc = BetaDiversityCalc(dist_unweighted_unifrac, 'unifrac', True)
        matrix, labels = beta_calc(data_path=self.l19_fp,
                                   tree_path=self.l19_tree, result_path=None, log_path=None)
        self.assertEqual(labels, self.l19_sample_names)

    def test_l19_unifrac_escaped_names(self):
        """beta calc works for unifrac when tips names are escaped in newick
        """
        beta_calc = BetaDiversityCalc(dist_unweighted_unifrac, 'unifrac', True)
        non_escaped_result = beta_calc(data_path=self.l19_fp,
                                       tree_path=self.l19_tree, result_path=None, log_path=None)

        l19_tree_str = "(((('tax7':0.1,'tax3':0.2):.98,tax8:.3, 'tax4':.3):.4,\
 (('ta_x1':0.3, tax6:.09):0.43,tax2:0.4):0.5):.2, (tax9:0.3, 'endbigtaxon':.08));"

        fd, tree_fp = mkstemp(prefix='Beta_div_tests', suffix='.tre')
        os.close(fd)
        open(tree_fp, 'w').write(l19_tree_str)
        self.files_to_remove.append(tree_fp)
        escaped_result = beta_calc(data_path=self.l19_w_underscore_fp,
                                   tree_path=tree_fp, result_path=None,
                                   log_path=None)
        npt.assert_almost_equal(escaped_result[0], non_escaped_result[0])
        self.assertEqual(escaped_result[1], non_escaped_result[1])

    def single_file_beta(
            self, otu_table_string, tree_string, missing_sams=None,
            use_metric_list=False):
        """ running single_file_beta should give same result using --rows"""
        if missing_sams is None:
            missing_sams = []
        # setup
        fd, input_path = mkstemp(suffix='.txt')
        os.close(fd)
        in_fname = os.path.split(input_path)[1]
        f = open(input_path, 'w')
        f.write(otu_table_string)
        f.close()
        fd, tree_path = mkstemp(suffix='.tre')
        os.close(fd)
        f = open(tree_path, 'w')
        f.write(tree_string)
        f.close()
        metrics = list_known_nonphylogenetic_metrics()
        metrics.extend(list_known_phylogenetic_metrics())
        output_dir = mkdtemp()

        # new metrics that don't trivially parallelize must be dealt with
        # carefully
        warnings.filterwarnings('ignore', 'dissimilarity binary_dist_chisq is\
 not parallelized, calculating the whole matrix...')
        warnings.filterwarnings('ignore', 'dissimilarity dist_chisq is not\
 parallelized, calculating the whole matrix...')
        warnings.filterwarnings('ignore', 'dissimilarity dist_gower is not\
 parallelized, calculating the whole matrix...')
        warnings.filterwarnings('ignore', 'dissimilarity dist_hellinger is\
 not parallelized, calculating the whole matrix...')
        warnings.filterwarnings('ignore', 'unifrac had no information for\
 sample M*')

        self.files_to_remove.extend([input_path, tree_path])
        self.folders_to_remove.append(output_dir)
        os.mkdir(output_dir + '/ft/')

        for metric in metrics:
            # do it
            if use_metric_list:
                single_file_beta(input_path, [metric], tree_path, output_dir,
                                 rowids=None)
            else:
                single_file_beta(input_path, metric, tree_path, output_dir,
                                 rowids=None)
            sams, dmtx = parse_distmat(open(output_dir + '/' +
                                            metric + '_' + in_fname))

            # do it by rows
            for i in range(len(sams)):
                if sams[i] in missing_sams:
                    continue
                rows = sams[i]
                row_outname = output_dir + '/' + metric + '_' +\
                    in_fname
                if use_metric_list:
                    single_file_beta(input_path, [metric], tree_path,
                                     output_dir, rowids=rows)
                else:
                    single_file_beta(input_path, metric, tree_path, output_dir,
                                     rowids=rows)
                col_sams, row_sams, row_dmtx = parse_matrix(open(row_outname))

                self.assertEqual(row_dmtx.shape, (len(rows.split(',')),
                                                  len(sams)))

                # make sure rows same as full
                for j in range(len(rows.split(','))):
                    for k in range(len(sams)):
                        row_v1 = row_dmtx[j, k]
                        full_v1 =\
                            dmtx[sams.index(row_sams[j]),
                                 sams.index(col_sams[k])]
                        npt.assert_almost_equal(row_v1, full_v1)

            # full tree run:
            if 'full_tree' in str(metric).lower():
                continue
            # do it by rows with full tree
            for i in range(len(sams)):
                if sams[i] in missing_sams:
                    continue
                rows = sams[i]

                row_outname = output_dir + '/ft/' + metric + '_' +\
                    in_fname
                if use_metric_list:
                    single_file_beta(input_path, [metric], tree_path,
                                     output_dir + '/ft/', rowids=rows, full_tree=True)
                else:
                    single_file_beta(input_path, metric, tree_path,
                                     output_dir + '/ft/', rowids=rows, full_tree=True)
                col_sams, row_sams, row_dmtx = parse_matrix(open(row_outname))

                self.assertEqual(row_dmtx.shape, (len(rows.split(',')),
                                                  len(sams)))

                # make sure rows same as full
                for j in range(len(rows.split(','))):
                    for k in range(len(sams)):
                        row_v1 = row_dmtx[j, k]
                        full_v1 =\
                            dmtx[sams.index(row_sams[j]),
                                 sams.index(col_sams[k])]
                        npt.assert_almost_equal(row_v1, full_v1)

            # do it with full tree
            if use_metric_list:
                single_file_beta(input_path, [metric], tree_path,
                                 output_dir + '/ft/', rowids=None, full_tree=True)
            else:
                single_file_beta(input_path, metric, tree_path,
                                 output_dir + '/ft/', rowids=None, full_tree=True)
            sams_ft, dmtx_ft = parse_distmat(open(output_dir + '/ft/' +
                                                  metric + '_' + in_fname))
            self.assertEqual(sams_ft, sams)
            npt.assert_almost_equal(dmtx_ft, dmtx)

    def test_single_file_beta(self):
        self.single_file_beta(l19_otu_table, l19_tree)

    def test_single_file_beta_metric_list(self):
        self.single_file_beta(l19_otu_table, l19_tree, use_metric_list=True)

    def test_single_file_beta_missing(self):
        self.single_file_beta(missing_otu_table, missing_tree,
                              missing_sams=['M'])

    def test_single_file_beta_missin_metric_list(self):
        self.single_file_beta(missing_otu_table, missing_tree,
                              missing_sams=['M'], use_metric_list=True)


l19_otu_table = """{"rows": [{"id": "tax1", "metadata": {}}, {"id": "tax2",\
 "metadata": {}}, {"id": "tax3", "metadata": {}}, {"id": "tax4", "metadata":\
 {}}, {"id": "endbigtaxon", "metadata": {}}, {"id": "tax6", "metadata": {}},\
 {"id": "tax7", "metadata": {}}, {"id": "tax8", "metadata": {}}, {"id": "tax9",\
 "metadata": {}}], "format": "Biological Observation Matrix v0.9", "data":\
 [[0, 0, 7.0], [0, 1, 4.0], [0, 2, 2.0], [0, 3, 1.0], [1, 0, 1.0], [1, 1, 2.0],\
 [1, 2, 4.0], [1, 3, 7.0], [1, 4, 8.0], [1, 5, 7.0], [1, 6, 4.0], [1, 7, 2.0],\
 [1, 8, 1.0], [2, 5, 1.0], [2, 6, 2.0], [2, 7, 4.0], [2, 8, 7.0], [2, 9, 8.0],\
 [2, 10, 7.0], [2, 11, 4.0], [2, 12, 2.0], [2, 13, 1.0], [3, 10, 1.0],\
 [3, 11, 2.0], [3, 12, 4.0], [3, 13, 7.0], [3, 14, 8.0], [3, 15, 7.0],\
 [3, 16, 4.0], [3, 17, 2.0], [3, 18, 1.0], [4, 15, 1.0], [4, 16, 2.0],\
 [4, 17, 4.0], [4, 18, 7.0], [5, 1, 1.0], [5, 2, 1.0], [6, 6, 2.0],\
 [6, 7, 1.0], [7, 11, 3.0], [7, 12, 1.0], [8, 16, 4.0], [8, 17, 1.0]],\
 "columns": [{"id": "sam1", "metadata": null}, {"id": "sam2", "metadata":\
 null}, {"id": "sam3", "metadata": null}, {"id": "sam4", "metadata": null},\
 {"id": "sam5", "metadata": null}, {"id": "sam6", "metadata": null}, {"id":\
 "sam7", "metadata": null}, {"id": "sam8", "metadata": null}, {"id": "sam9",\
 "metadata": null}, {"id": "sam_middle", "metadata": null}, {"id": "sam11",\
 "metadata": null}, {"id": "sam12", "metadata": null}, {"id": "sam13",\
 "metadata": null}, {"id": "sam14", "metadata": null}, {"id": "sam15",\
 "metadata": null}, {"id": "sam16", "metadata": null}, {"id": "sam17",\
 "metadata": null}, {"id": "sam18", "metadata": null}, {"id": "sam19",\
 "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2520",\
 "matrix_type": "sparse", "shape": [9, 19], "format_url":\
 "http://www.qiime.org/svn_documentation/documentation/biom_format.html",\
 "date": "2011-12-20T19:03:28.130403", "type": "OTU table", "id": null,\
 "matrix_element_type": "float"}"""

l19_tree = """((((tax7:0.1,tax3:0.2):.98,tax8:.3, tax4:.3):.4, ((tax1:0.3,\
 tax6:.09):0.43,tax2:0.4):0.5):.2, (tax9:0.3, endbigtaxon:.08));"""

missing_tree = """((a:1,b:2):4,((c:3, (j:1,k:2)mt:17),(d:1,e:1):2):3)"""
missing_otu_table = """{"rows": [{"id": "a", "metadata": {}}, {"id": "c",\
 "metadata": {}}, {"id": "b", "metadata": {}}, {"id": "e", "metadata": {}},\
 {"id": "m", "metadata": {}}, {"id": "d", "metadata": {}}], "format":\
 "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 2, 2.0],\
 [1, 1, 1.0], [2, 0, 1.0], [2, 1, 1.0], [3, 2, 1.0], [4, 3, 88.0],\
 [5, 1, 3.0]], "columns": [{"id": "A", "metadata": null}, {"id": "B",\
 "metadata": null}, {"id": "C", "metadata": null}, {"id": "M", "metadata":\
 null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2520", "matrix_type":\
 "sparse", "shape": [6, 4], "format_url":\
 "http://www.qiime.org/svn_documentation/documentation/biom_format.html",\
 "date": "2011-12-20T19:04:03.875636", "type": "OTU table", "id": null,\
 "matrix_element_type": "float"}"""

# run tests if called from command line
if __name__ == '__main__':
    main()
