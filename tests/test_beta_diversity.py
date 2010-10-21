#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"


"""Contains tests for performing beta diversity analyses."""

import numpy
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from cogent.app.util import get_tmp_filename
from cogent.core.tree import PhyloNode
from cogent.maths.distance_transform import dist_chisq
from qiime.parse import parse_newick, parse_distmat, parse_matrix
from qiime.format import format_otu_table
from qiime.beta_diversity import BetaDiversityCalc, single_file_beta,\
list_known_nonphylogenetic_metrics, list_known_phylogenetic_metrics
from qiime.beta_metrics import dist_unweighted_unifrac

import os
import shutil

class BetaDiversityCalcTests(TestCase):
    """Tests of the BetaDiversityCalc class"""
    def setUp(self):
        self.l19_data = numpy.array([
            [7,1,0,0,0,0,0,0,0],
            [4,2,0,0,0,1,0,0,0],
            [2,4,0,0,0,1,0,0,0],
            [1,7,0,0,0,0,0,0,0],
            [0,8,0,0,0,0,0,0,0],
            [0,7,1,0,0,0,0,0,0],
            [0,4,2,0,0,0,2,0,0],
            [0,2,4,0,0,0,1,0,0],
            [0,1,7,0,0,0,0,0,0],
            [0,0,8,0,0,0,0,0,0],
            [0,0,7,1,0,0,0,0,0],
            [0,0,4,2,0,0,0,3,0],
            [0,0,2,4,0,0,0,1,0],
            [0,0,1,7,0,0,0,0,0],
            [0,0,0,8,0,0,0,0,0],
            [0,0,0,7,1,0,0,0,0],
            [0,0,0,4,2,0,0,0,4],
            [0,0,0,2,4,0,0,0,1],
            [0,0,0,1,7,0,0,0,0]
            ])
        self.l19_sample_names = ['sam1', 'sam2', 'sam3', 'sam4', 'sam5','sam6',\
        'sam7', 'sam8', 'sam9', 'sam_middle', 'sam11', 'sam12', 'sam13', \
        'sam14', 'sam15', 'sam16', 'sam17', 'sam18', 'sam19']
        self.l19_taxon_names =  ['tax1', 'tax2', 'tax3', 'tax4', 'endbigtaxon',\
        'tax6', 'tax7', 'tax8', 'tax9']
        self.l19_taxon_names_w_underscore =  ['ta_x1', 'tax2', 'tax3', 'tax4', 
         'endbigtaxon', 'tax6', 'tax7', 'tax8', 'tax9']

        self.l19_str = format_otu_table(
            self.l19_sample_names, 
            self.l19_taxon_names, 
            self.l19_data.T)
            
        self.l19_str_w_underscore = format_otu_table(
            self.l19_sample_names, 
            self.l19_taxon_names_w_underscore, 
            self.l19_data.T)

        self.l19_tree_str = '((((tax7:0.1,tax3:0.2):.98,tax8:.3, tax4:.3):.4, ((tax1:0.3, tax6:.09):0.43,tax2:0.4):0.5):.2, (tax9:0.3, endbigtaxon:.08));'
        self.l19_tree = parse_newick(self.l19_tree_str, PhyloNode)
        
        self.files_to_remove = []
        self.folders_to_remove = []
        
    def tearDown(self):
        remove_files(self.files_to_remove)
        for folder in self.folders_to_remove:
            shutil.rmtree(folder)

    def test_l19_chi(self):
        """beta calc run should return same values as directly calling metric"""
        beta_calc_chisq = BetaDiversityCalc(dist_chisq, 'chi square', False)
        matrix, labels = beta_calc_chisq(data_path=self.l19_str, tree_path=None)
        self.assertEqual(labels, self.l19_sample_names)
        self.assertFloatEqual(matrix, dist_chisq(self.l19_data))
        
    
    def test_l19_unifrac(self):
        """beta calc run should also work for phylo metric"""
        beta_calc = BetaDiversityCalc(dist_unweighted_unifrac, 'unifrac', True)
        data_path = self.l19_str
        matrix, labels = beta_calc(data_path=data_path, 
            tree_path=self.l19_tree, result_path=None, log_path=None)
        self.assertEqual(labels, self.l19_sample_names)
        
    def test_l19_unifrac_escaped_names(self):
        """beta calc works for unifrac when tips names are escaped in newick 
        """
        beta_calc = BetaDiversityCalc(dist_unweighted_unifrac, 'unifrac', True)
        non_escaped_result = beta_calc(data_path=self.l19_str, 
            tree_path=self.l19_tree, result_path=None, log_path=None)
            
        l19_tree_str = "(((('tax7':0.1,'tax3':0.2):.98,tax8:.3, 'tax4':.3):.4, (('ta_x1':0.3, tax6:.09):0.43,tax2:0.4):0.5):.2, (tax9:0.3, 'endbigtaxon':.08));"
        
        tree_fp = get_tmp_filename(prefix='Beta_div_tests',suffix='.tre')
        open(tree_fp,'w').write(l19_tree_str)
        self.files_to_remove.append(tree_fp)
        escaped_result = beta_calc(data_path=self.l19_str_w_underscore, 
            tree_path=tree_fp, result_path=None, log_path=None)
        self.assertEqual(escaped_result,non_escaped_result)
        
    def test_single_file_beta(self):
        """ running single_file_beta should give same result using --rows"""
        # setup
        input_path = get_tmp_filename()
        in_fname = os.path.split(input_path)[1]
        f = open(input_path,'w')
        f.write(l19_otu_table)
        f.close()
        tree_path = get_tmp_filename()
        f = open(tree_path,'w')
        f.write(l19_tree)
        f.close()
        metrics = list_known_nonphylogenetic_metrics()
        metrics.extend(list_known_phylogenetic_metrics())
        output_dir = get_tmp_filename(suffix = '')
        os.mkdir(output_dir)        
        
        self.files_to_remove.extend([input_path,tree_path])
        self.folders_to_remove.append(output_dir)
        
        for metric in metrics:
            # do it
            single_file_beta(input_path, ','.join(metrics), tree_path, output_dir,
                rowids=None)
            sams, dmtx = parse_distmat(open(output_dir + '/' +\
                metric + '_' + in_fname))

            # do it by rows
            rows = 'sam5,sam3'
            row_outname = output_dir + '/' + metric + '_' +\
                '_'.join(rows.split(',')) + '_' + in_fname
            single_file_beta(input_path, ','.join(metrics), tree_path, output_dir,
                rowids=rows)
            col_sams, row_sams, row_dmtx = parse_matrix(open(row_outname))
            
            self.assertEqual(row_dmtx.shape, (len(rows.split(',')),len(sams)))
        
            # make sure rows same as full
            for i in range(len(rows.split(','))):
                for j in range(len(sams)):
                    row_v1 = row_dmtx[i,j]
                    full_v1 =\
                        dmtx[sams.index(row_sams[i]),sams.index(col_sams[j])]
                    self.assertFloatEqual(row_v1, full_v1)
        

l19_otu_table = """#comment line
#OTU ID	sam1	sam2	sam3	sam4	sam5	sam6	sam7	sam8	sam9	sam_middle	sam11	sam12	sam13	sam14	sam15	sam16	sam17	sam18	sam19
tax1	7	4	2	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
tax2	1	2	4	7	8	7	4	2	1	0	0	0	0	0	0	0	0	0	0
tax3	0	0	0	0	0	1	2	4	7	8	7	4	2	1	0	0	0	0	0
tax4	0	0	0	0	0	0	0	0	0	0	1	2	4	7	8	7	4	2	1
endbigtaxon	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	2	4	7
tax6	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
tax7	0	0	0	0	0	0	2	1	0	0	0	0	0	0	0	0	0	0	0
tax8	0	0	0	0	0	0	0	0	0	0	0	3	1	0	0	0	0	0	0
tax9	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	4	1	0"""

l19_tree = """((((tax7:0.1,tax3:0.2):.98,tax8:.3, tax4:.3):.4, ((tax1:0.3, tax6:.09):0.43,tax2:0.4):0.5):.2, (tax9:0.3, endbigtaxon:.08));"""


#run tests if called from command line
if __name__ == '__main__':
    main()
