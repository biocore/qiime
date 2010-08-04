#!/usr/bin/env python

"""Tests of code for performing significance tests of OTU/category associations
"""

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Catherine Lozupone", "Dan Knights"] 
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from qiime.otu_category_significance import filter_OTUs, \
    make_contingency_matrix, run_single_G_test, run_G_test_OTUs, \
    add_fdr_correction_to_results, output_results_G_test, \
    run_single_ANOVA, run_ANOVA_OTUs, output_results_ANOVA,\
    run_correlation_OTUs, run_single_correlation, output_results_correlation,\
    parse_otu_table, parse_category_mapping,\
    aggregate_multiple_results_ANOVA, aggregate_multiple_results_G_test,\
    aggregate_multiple_results_correlation, get_common_OTUs,\
    test_wrapper_multiple, test_wrapper
from numpy import array
from qiime.otu_category_significance import parse_otu_table
from cogent.util.dict2d import Dict2D
from cogent.app.util import get_tmp_filename
from os import remove

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_filter_OTUs(self):
        """filter_OTUs works"""
        otu_table = """#Full OTU Counts
#OTU ID\tsample1\tsample2\tsample3
0\t0\t2\t0
1\t1\t0\t0
2\t1\t1\t1""".split('\n')
        OTU_sample_info, num_samples, taxonomy_info = parse_otu_table(otu_table)
        result = filter_OTUs(OTU_sample_info, 2, 3)
        self.assertEqual(result, [])
        result = filter_OTUs(OTU_sample_info, 1, 3)
        self.assertEqual(result, ['1', '0'])
        
        result = filter_OTUs(OTU_sample_info, 2, 3, False)
        self.assertEqual(result, ['2'])
        result = filter_OTUs(OTU_sample_info, 1, 3, False)
        self.assertEqual(result, ['1', '0', '2'])
        #test that is works if a category mapping file is supplied
        cat_mapping = {'sample2': '0', 'sample3': '1'}
        result = filter_OTUs(OTU_sample_info, 1,\
                        category_mapping_info=cat_mapping)
        self.assertEqual(result, ['0'])

    def test_make_contingency_matrix(self):
        """make_contingency_matrix works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}
        OTU_sample_info = {'0': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '1'},
        '1': {'sample1': '1', 'sample2': '1', 'sample3': '0', 'sample4': '0'},
        '2': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '0'},
        '3': {'sample1': '0', 'sample2': '0', 'sample3': '0', 'sample4': '1'}}
        category_values = ['A', 'B', 'C']
        result = make_contingency_matrix('0', category_info, OTU_sample_info, category_values)
        self.assertEqual(result['OTU_pos']['B_pos'], 1)
        self.assertEqual(result['OTU_pos']['C_pos'], 1)
        self.assertEqual(result['OTU_pos']['A_pos'], 0)
        self.assertEqual(result['OTU_neg']['A_pos'], 2)
        self.assertEqual(result['OTU_neg']['B_pos'], 0)
        self.assertEqual(result['OTU_neg']['C_pos'], 0)

    def test_run_single_G_test(self):
        """run_single_G_test works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}
        OTU_sample_info = {'0': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '1'},
        '1': {'sample1': '1', 'sample2': '1', 'sample3': '0', 'sample4': '0'},
        '2': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '0'},
        '3': {'sample1': '0', 'sample2': '0', 'sample3': '0', 'sample4': '1'}}
        category_values = ['A', 'B', 'C']
        g_val, prob, contingency_matrix = run_single_G_test('0', category_info, OTU_sample_info, category_values)
        self.assertFloatEqual(g_val, 4.29304060218)
        self.assertFloatEqual(prob, 0.508041627088)
        self.assertEqual(contingency_matrix, {'OTU_pos': {'B_pos': [1, 0.5], 'C_pos': [1, 0.5], 'A_pos': [0, 1.0]}, 'OTU_neg': {'B_pos': [0, 0.5], 'C_pos': [0, 0.5], 'A_pos': [2, 1.0]}})


    def test_run_single_ANOVA(self):
        """run_single_ANOVA works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'B'}
        OTU_sample_info = {'0': {'sample1': '5', 'sample2': '10', 'sample3': '2', 'sample4': '1'},
        '1': {'sample1': '1', 'sample2': '1', 'sample3': '0', 'sample4': '0'},
        '2': {'sample1': '2', 'sample2': '0', 'sample3': '1', 'sample4': '0'},
        '3': {'sample1': '0', 'sample2': '0', 'sample3': '0', 'sample4': '1'}}
        category_values = ['A', 'B']
        group_means, prob = run_single_ANOVA('0', category_info,\
            OTU_sample_info, category_values)
        self.assertEqual(group_means, [7.5, 1.5])
        self.assertFloatEqual(prob, 0.142857142857)

    def test_run_single_correlation(self):
        """run_single_correlation works"""
        category_info = {'sample1': '1',
                        'sample2': '2',
                        'sample3': '3',
                        'sample4': '4'}
        OTU_sample_info = {'0': {'sample1': '1', 'sample2': '2', 'sample3': '2', 'sample4': '4'},
        '1': {'sample1': '1', 'sample2': '2', 'sample3': '3', 'sample4': '4'},
        '2': {'sample1': '4', 'sample2': '3', 'sample3': '2', 'sample4': '1'},
        '3': {'sample1': '1', 'sample2': '2', 'sample3': '3', 'sample4': '5'}}
        category_values = ['A', 'B']
        r, prob = run_single_correlation('0', category_info,\
            OTU_sample_info)
        self.assertFloatEqual(r, 0.923380516877)
        self.assertFloatEqual(prob, 0.0766194831234)

    def test_run_ANOVA_OTUs(self):
        """run_ANOVA_OTUs works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'B'}
        OTU_sample_info = {'0': {'sample1': '5', 'sample2': '10', 'sample3': '2', 'sample4': '1'},
        '1': {'sample1': '1', 'sample2': '0', 'sample3': '0', 'sample4': '2'},
        '2': {'sample1': '2', 'sample2': '1', 'sample3': '10', 'sample4': '15'},
        '3': {'sample1': '1', 'sample2': '1.5', 'sample3': '1.4', 'sample4': '1.3'}}
        category_values = ['A', 'B']
        result = run_ANOVA_OTUs(['0', '1', '3'], category_info,\
            OTU_sample_info, category_values)
        self.assertEqual(result, {'1': [[0.5, 1.0], 0.69848865542223582, 2.0954659662667074], '0': [[7.5, 1.5], 0.14285714285714285, 0.42857142857142855], '3': [[1.25, 1.3500000000000001], 0.73273875808757438, 2.1982162742627231]})

    def test_output_results_ANOVA(self):
        """output_results_ANOVA works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'B'}
        OTU_sample_info = {'0': {'sample1': '5', 'sample2': '10', 'sample3': '2', 'sample4': '1'},
        '1': {'sample1': '1', 'sample2': '0', 'sample3': '0', 'sample4': '2'},
        '2': {'sample1': '2', 'sample2': '1', 'sample3': '10', 'sample4': '15'},
        '3': {'sample1': '1', 'sample2': '1.5', 'sample3': '1.4', 'sample4': '1.3'},
        '4': {'sample1': '20', 'sample2': '16', 'sample3': '1.4', 'sample4': '1.3'}}
        category_values = ['A', 'B']
        taxonomy_info = {'0': 'taxon1',
                        '1': 'taxon2',
                        '2': 'taxon3',
                        '3': 'taxon4',
                        '4': 'taxon5'}
        ANOVA_results = run_ANOVA_OTUs(['0', '2', '1', '3', '4'], category_info,\
            OTU_sample_info, category_values)
        output = output_results_ANOVA(ANOVA_results, category_values,\
                        taxonomy_info)
        self.assertEqual(output, ['OTU\tprob\tBonferroni_corrected\tFDR_corrected\tA_mean\tB_mean\tConsensus Lineage', '1\t0.698488655422\t3.49244327711\t0.873110819278\t0.5\t1.0\ttaxon2', '0\t0.142857142857\t0.714285714286\t0.238095238095\t7.5\t1.5\ttaxon1', '3\t0.732738758088\t3.66369379044\t0.732738758088\t1.25\t1.35\ttaxon4', '2\t0.0497447318605\t0.248723659303\t0.124361829651\t1.5\t12.5\ttaxon3', '4\t0.0141325222337\t0.0706626111683\t0.0706626111683\t18.0\t1.35\ttaxon5'])
    
    def test_run_G_test_OTUs(self):
        """run_G_test_OTUs works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}
        OTU_sample_info = {'0': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '1'},
        '1': {'sample1': '1', 'sample2': '1', 'sample3': '0', 'sample4': '0'},
        '2': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '0'},
        '3': {'sample1': '0', 'sample2': '0', 'sample3': '0', 'sample4': '1'}}
        category_values = ['A', 'B', 'C']
        result = run_G_test_OTUs(['0', '1', '3'], category_info, \
            OTU_sample_info, category_values)
        self.assertEqual(result.keys(), ['1', '0', '3'])
        self.assertFloatEqual(result['0'][3], 1.5241248812628101)
        self.assertFloatEqual(result['0'][0], 4.29304060218)
        self.assertFloatEqual(result['0'][1], 0.508041627088)
        self.assertEqual(result['0'][2], {'OTU_pos': {'B_pos': [1, 0.5], 'C_pos': [1, 0.5], 'A_pos': [0, 1.0]}, 'OTU_neg': {'B_pos': [0, 0.5], 'C_pos': [0, 0.5], 'A_pos': [2, 1.0]}})

    def test_fdr_correction_G_test(self):
        """fdr_correction_G_test works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}
        OTU_sample_info = {'0': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '1'},
        '1': {'sample1': '1', 'sample2': '1', 'sample3': '1', 'sample4': '0'},
        '2': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '1'},
        '3': {'sample1': '0', 'sample2': '1', 'sample3': '0', 'sample4': '1'}}
        category_values = ['A', 'B', 'C']
        G_test_results = run_G_test_OTUs(['0', '1', '3'], category_info, \
            OTU_sample_info, category_values)
        G_test_results = add_fdr_correction_to_results(G_test_results)
        self.assertFloatEqual(G_test_results['0'][-1], 1.52412488126)
        self.assertFloatEqual(G_test_results['1'][-1], 0.938976340277)
        self.assertFloatEqual(G_test_results['3'][-1], 0.828522198394)

    def test_output_results_G_test(self):
        """output_results works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}
        OTU_sample_info = {'0': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '1'},
        '1': {'sample1': '1', 'sample2': '1', 'sample3': '1', 'sample4': '0'},
        '2': {'sample1': '0', 'sample2': '0', 'sample3': '1', 'sample4': '1'},
        '3': {'sample1': '0', 'sample2': '1', 'sample3': '0', 'sample4': '1'}}
        category_values = ['A', 'B', 'C']
        taxonomy_info = {'0': 'taxon1',
                        '1': 'taxon2',
                        '2': 'taxon3',
                        '3': 'taxon4'}
        G_test_results = run_G_test_OTUs(['0', '1', '3'], category_info, \
            OTU_sample_info, category_values)
        output = output_results_G_test(G_test_results, {})
        self.assertEqual(output, ['OTU\tg_val\tg_prob\tBonferroni_corrected\tFDR_corrected\tOTU_pos##B_pos\tOTU_pos##C_pos\tOTU_pos##A_pos\tOTU_neg##B_pos\tOTU_neg##C_pos\tOTU_neg##A_pos', '1\t3.48284992796\t0.625984226851\t1.87795268055\t0.938976340277\t[1, 0.75]\t[0, 0.75]\t[2, 1.5]\t[0, 0.25]\t[1, 0.25]\t[0, 0.5]', '0\t4.29304060218\t0.508041627088\t1.52412488126\t1.52412488126\t[1, 0.5]\t[1, 0.5]\t[0, 1.0]\t[0, 0.5]\t[0, 0.5]\t[2, 1.0]', '3\t2.14652030109\t0.828522198394\t2.48556659518\t0.828522198394\t[0, 0.5]\t[1, 0.5]\t[1, 1.0]\t[1, 0.5]\t[0, 0.5]\t[1, 1.0]'])
        output = output_results_G_test(G_test_results, taxonomy_info)
        self.assertEqual(output, ['OTU\tg_val\tg_prob\tBonferroni_corrected\tFDR_corrected\tOTU_pos##B_pos\tOTU_pos##C_pos\tOTU_pos##A_pos\tOTU_neg##B_pos\tOTU_neg##C_pos\tOTU_neg##A_pos\tConsensus Lineage', '1\t3.48284992796\t0.625984226851\t1.87795268055\t0.938976340277\t[1, 0.75]\t[0, 0.75]\t[2, 1.5]\t[0, 0.25]\t[1, 0.25]\t[0, 0.5]\ttaxon2', '0\t4.29304060218\t0.508041627088\t1.52412488126\t1.52412488126\t[1, 0.5]\t[1, 0.5]\t[0, 1.0]\t[0, 0.5]\t[0, 0.5]\t[2, 1.0]\ttaxon1', '3\t2.14652030109\t0.828522198394\t2.48556659518\t0.828522198394\t[0, 0.5]\t[1, 0.5]\t[1, 1.0]\t[1, 0.5]\t[0, 0.5]\t[1, 1.0]\ttaxon4'])

    def test_run_correlation_OTUs(self):
        """run_single_correlation works"""
        category_info = {'sample1':'0.1',
                        'sample2':'0.2',
                        'sample3':'0.3',
                        'sample4':'0.4'}
        OTU_sample_info = {'0': {'sample1': '0', 'sample2': '1', 'sample3': '2', 'sample4': '3'},
        '1': {'sample1': '7', 'sample2': '5', 'sample3': '3', 'sample4': '1'},
        '2': {'sample1': '4', 'sample2': '4.2', 'sample3': '4', 'sample4': '4'},
        '3': {'sample1': '0', 'sample2': '1', 'sample3': '0', 'sample4': '1'}}
        OTU_list = ['1', '0', '2', '3'] 
        result = run_correlation_OTUs(OTU_list, category_info, OTU_sample_info)
        #OTU 0 should be positively correlated, 1 negative, and 2&3 neutral
        self.assertFloatEqual(result['0'], [1.0, 0.0, 0.0])
        self.assertFloatEqual(result['1'], [-0.99999999999999956, 4.4408920985006281e-16, 1.7763568394002513e-15])
        self.assertFloatEqual(result['2'], [-0.25819888974715061, 0.74180111025284945, 2.9672044410113978])
        self.assertFloatEqual(result['3'], [0.44721359549995815, 0.55278640450004191, 2.2111456180001676])

        #test that appropriate error is raised is categorical
        category_info = {'sample1':'A',
                        'sample2':'B',
                        'sample3':'A',
                        'sample4':'B'}
        self.assertRaises(ValueError, run_correlation_OTUs, OTU_list, category_info, OTU_sample_info)
        
    def test_output_results_correlation(self):
        """output_results_correlation works"""
        category_info = {'sample1':'0.1',
                        'sample2':'0.2',
                        'sample3':'0.3',
                        'sample4':'0.4'}
        OTU_sample_info = {'0': {'sample1': '0', 'sample2': '1', 'sample3': '2', 'sample4': '3'},
        '1': {'sample1': '7', 'sample2': '5', 'sample3': '3', 'sample4': '1'},
        '2': {'sample1': '4', 'sample2': '4.2', 'sample3': '4', 'sample4': '4'},
        '3': {'sample1': '0', 'sample2': '1', 'sample3': '0', 'sample4': '1'}}
        taxonomy_info = {'0': 'taxon1',
                        '1': 'taxon2',
                        '2': 'taxon3',
                        '3': 'taxon4'}
        OTU_list = ['1', '0', '2', '3'] 
        result = run_correlation_OTUs(OTU_list, category_info, OTU_sample_info)
        output = output_results_correlation(result, taxonomy_info)
        self.assertEqual(output[0], 'OTU\tprob\tBonferroni_corrected\tFDR_corrected\tr\tConsensus Lineage')
        self.assertEqual(output[1], '1\t4.4408920985e-16\t1.7763568394e-15\t8.881784197e-16\t-1.0\ttaxon2')
        self.assertEqual(len(output), 5)

    def test_parse_otu_table(self):
        """parse otu_table works"""
        otu_table = """#Full OTU Counts
#OTU ID\tsample1\tsample2\tsample3
0\t0\t2\t0
1\t1\t0\t0
2\t1\t1\t1""".split('\n')
        result, num_samples, taxonomy_info = parse_otu_table(otu_table)
        self.assertEqual(result['1'], {'sample1': '1', 'sample3': '0', 'sample2': '0'})
        self.assertEqual(result['0'], {'sample1': '0', 'sample3': '0', 'sample2': '2'})
        self.assertEqual(result['2'], {'sample1': '1', 'sample3': '1', 'sample2': '1'})
        self.assertEqual(num_samples, 3)
        self.assertEqual(taxonomy_info, {})

        #test that it parses otu tables with taxonomy fields appropriately
        otu_table = """#Full OTU Counts
#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage
0\t0\t2\t0\tBacteria; Bacteroidetes; Bacteroidales; Parabacteroidaceae; Unclassified; otu_475
1\t1\t0\t0\tBacteria; Bacteroidetes; Bacteroidales; adhufec77-25; Barnesiella; Barnesiella_viscericola; otu_369
2\t1\t1\t1\tBacteria; Firmicutes; Clostridia; Clostridiales; Faecalibacterium; Unclassified; otu_1121""".split('\n')
        result, num_samples, taxonomy_info = parse_otu_table(otu_table)
        self.assertEqual(result['1'], {'sample1': '1', 'sample3': '0', 'sample2': '0'})
        self.assertEqual(result['0'], {'sample1': '0', 'sample3': '0', 'sample2': '2'})
        self.assertEqual(result['2'], {'sample1': '1', 'sample3': '1', 'sample2': '1'})
        self.assertEqual(num_samples, 3)
        self.assertEqual(taxonomy_info, {'1': 'Bacteria; Bacteroidetes; Bacteroidales; adhufec77-25; Barnesiella; Barnesiella_viscericola; otu_369', '0': 'Bacteria; Bacteroidetes; Bacteroidales; Parabacteroidaceae; Unclassified; otu_475', '2': 'Bacteria; Firmicutes; Clostridia; Clostridiales; Faecalibacterium; Unclassified; otu_1121'})

    def test_parse_category_mapping(self):
        """parse_category_mapping works"""
        category_mapping = """#SampleID\tcat1\tcat2
sample1\tA\t0
sample2\tB\t8.0
sample3\tC\t1.0""".split('\n')
        result, cat_vals = parse_category_mapping(category_mapping, 'cat1')
        self.assertEqual(result, {'sample1': 'A', 'sample3': 'C', 'sample2': 'B'})
        self.assertEqual(cat_vals, (['A', 'B', 'C']))
        result, cat_vals = parse_category_mapping(category_mapping, 'cat2', threshold=5.0)
        self.assertEqual(result, {'sample1': '0', 'sample3': '0', 'sample2': '1'})
        self.assertEqual(cat_vals, (['0', '1']))
        
    def test_test_wrapper(self):
        """runs the specified statistical test"""
        #correlation

    def test_aggregate_multiple_results_ANOVA(self):
        """aggregate_multiple_results_ANOVA works"""
        all_results = {
            '0':[[[10,20],.05,.1], [[12,18],.1,.2]],
            '1':[[[15,14],.5,1], [[12,11],.6, 1.2]]
            }
        result = aggregate_multiple_results_ANOVA(all_results)
        self.assertFloatEqual(result['0'], [array([11,19]),.075, .15])
        self.assertFloatEqual(result['1'], [array([13.5,12.5]),.55, 1.1])
                       
    def test_aggregate_multiple_results_G_test(self):
        """aggregate_multiple_results_G_test works"""
        cm_0_0 = Dict2D({'OTU_pos':{'A_pos':2, 'B_pos':0},
                  'OTU_neg':{'A_pos':3,'B_pos':4}},Default=0, Pad=True)
        cm_0_1 = Dict2D({'OTU_pos':{'A_pos':5, 'B_pos':1},
                  'OTU_neg':{'A_pos':6,'B_pos':10}},Default=0, Pad=True)
        cm_1_0 = Dict2D({'OTU_pos':{'A_pos':3, 'B_pos':5},
                  'OTU_neg':{'A_pos':1,'B_pos':1}},Default=0, Pad=True)
        cm_1_1 = Dict2D({'OTU_pos':{'A_pos':4, 'B_pos':4},
                  'OTU_neg':{'A_pos':0,'B_pos':1}},Default=0, Pad=True)

        exp_cm_0 = Dict2D({'OTU_pos':{'A_pos':3.5, 'B_pos':0.5},
                  'OTU_neg':{'A_pos':4.5,'B_pos':7}},Default=0, Pad=True)
        exp_cm_1 = Dict2D({'OTU_pos':{'A_pos':3.5, 'B_pos':4.5},
                  'OTU_neg':{'A_pos':.5,'B_pos':1}},Default=0, Pad=True)
        all_results = {\
                 '0':[[.5, .05, cm_0_0, .1], [.6,.1,cm_0_1,.2]],
                 '1':[[.6, .5, cm_1_0, 1], [.8,.6, cm_1_1, 1.2]]
                 }
        results = aggregate_multiple_results_G_test(all_results)

        self.assertFloatEqual(results['0'], [0.55, 0.075, exp_cm_0, .15])
        self.assertFloatEqual(results['1'], [0.7, 0.55, exp_cm_1, 1.1])
        

    def test_aggregate_multiple_results_correlation(self):
        """aggregate_multiple_results_correlation works"""
        all_results = {
            '0':[[-1,.05,.1], [-.5,.1,.2]],
            '1':[[.4,.5,1], [.5,.6, 1.2]]
            }
        result = aggregate_multiple_results_correlation(all_results)
        self.assertFloatEqual(result['0'], [-.75, .075, .15])
        self.assertFloatEqual(result['1'], [.45,.55, 1.1])


    def test_get_common_OTUs(self):
        """get_common_OTUs works"""

        # create the temporary OTU tables
        otu_table1 = '\n'.join(['#Full OTU Counts',
                     '#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage',
                     '0\t0\t2\t0\tlineage0',
                     '1\t1\t0\t0\tlineage1',
                     '2\t1\t1\t1\tlineage2'])
        otu_table2 = '\n'.join(['#Full OTU Counts',
                     '#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage',
                     '0\t0\t2\t0\tlineage0',
                     '1\t1\t0\t0\tlineage1',
                     '2\t0\t1\t1\tlineage2'])
        otu_table3 = '\n'.join(['#Full OTU Counts',
                     '#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage',
                     '0\t0\t2\t0\tlineage0',
                     '2\t1\t1\t1\tlineage2'])
        category_info = {'sample1':'0.1',
                        'sample2':'0.2',
                        'sample3':'0.3'}
        OTU_list = ['1', '0', '2'] 

        fp1 = get_tmp_filename()
        fp2 = get_tmp_filename()
        fp3 = get_tmp_filename()
        try:
            f1 = open(fp1,'w')
            f2 = open(fp2,'w')
            f3 = open(fp3,'w')
        except IOError, e:
            raise e,"Could not create temporary files: %s, %s" %(f1,f2, f3)
        
        f1.write(otu_table1)
        f1.close()
        f2.write(otu_table2)
        f2.close()
        f3.write(otu_table3)
        f3.close()

        # case where one OTU is missing from one file
        otu_table_paths = [fp1,fp2,fp3]
        _filter = 0
        filter_all_samples = False
        otu_include = None
        OTU_list, taxonomy = get_common_OTUs(otu_table_paths, _filter, \
                                             category_info, \
                                             filter_all_samples, \
                                             otu_include=otu_include)
        exp_OTU_list = sorted(['0','2'])
        exp_taxonomy = {'0': 'lineage0', '2': 'lineage2'}
        self.assertEqual(sorted(OTU_list), exp_OTU_list)
        self.assertEqual(taxonomy,exp_taxonomy)

        # case where no OTUs should be filtered
        otu_table_paths = [fp1,fp2]
        _filter = 0
        filter_all_samples = False
        otu_include = None
        OTU_list, taxonomy = get_common_OTUs(otu_table_paths, _filter, \
                                             category_info, \
                                             filter_all_samples, \
                                             otu_include=otu_include)
        exp_OTU_list = sorted(['0','1','2'])
        exp_taxonomy = {'0': 'lineage0', '1':'lineage1', '2': 'lineage2'}
        self.assertEqual(sorted(OTU_list), exp_OTU_list)
        self.assertEqual(taxonomy,exp_taxonomy)

        # case where one OTU should be filtered due to filter_all_samples
        otu_table_paths = [fp1,fp2]
        _filter = 0
        filter_all_samples = True
        otu_include = None
        OTU_list, taxonomy = get_common_OTUs(otu_table_paths, _filter, \
                                             category_info, \
                                             filter_all_samples, \
                                             otu_include=otu_include)
        exp_OTU_list = sorted(['0','1'])
        exp_taxonomy = {'0': 'lineage0', '1':'lineage1'}
        self.assertEqual(sorted(OTU_list), exp_OTU_list)
        self.assertEqual(taxonomy,exp_taxonomy)

        # case where two OTUs should be filtered due to _filter value
        otu_table_paths = [fp1,fp2]
        _filter = 2
        filter_all_samples = False
        otu_include = None
        OTU_list, taxonomy = get_common_OTUs(otu_table_paths, _filter, \
                                             category_info, \
                                             filter_all_samples, \
                                             otu_include=otu_include)
        exp_OTU_list = sorted(['2'])
        exp_taxonomy = {'2':'lineage2'}
        self.assertEqual(sorted(OTU_list), exp_OTU_list)
        self.assertEqual(taxonomy,exp_taxonomy)


        # clean up temporary files
        remove(fp1)
        remove(fp2)
        remove(fp3)

    def test_test_wrapper_multiple(self):
        """test_wrapper_multiple works"""
        # create the temporary OTU tables
        otu_table1 = '\n'.join(['#Full OTU Counts',
                     '#OTU ID\tsample1\tsample2\tsample3\tsample4\tConsensus Lineage',
                     '0\t1\t2\t0\t1\tlineage0',
                     '1\t1\t0\t0\t1\tlineage1',
                     '2\t1\t1\t1\t1\tlineage2'])
        otu_table2 = '\n'.join(['#Full OTU Counts',
                     '#OTU ID\tsample1\tsample2\tsample3\tsample4\tConsensus Lineage',
                     '0\t0\t2\t0\t1\tlineage0',
                     '1\t1\t0\t0\t1\tlineage1',
                     '2\t0\t1\t1\t1\tlineage2'])
        otu_table3 = '\n'.join(['#Full OTU Counts',
                     '#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage',
                     '0\t0\t2\t0\t1\tlineage0',
                     '2\t1\t1\t1\t1\tlineage2'])
        category_mapping = ['#SampleID\tcat1\tcat2',
                                      'sample1\tA\t0',
                                      'sample2\tA\t8.0',
                                      'sample3\tB\t1.0',
                                      'sample4\tB\t1.0']
        OTU_list = ['1', '0'] 

        fp1 = get_tmp_filename()
        fp2 = get_tmp_filename()
        fp3 = get_tmp_filename()
        try:
            f1 = open(fp1,'w')
            f2 = open(fp2,'w')
            f3 = open(fp3,'w')
        except IOError, e:
            raise e,"Could not create temporary files: %s, %s" %(f1,f2, f3)
        
        f1.write(otu_table1)
        f1.close()
        f2.write(otu_table2)
        f2.close()
        f3.write(otu_table3)
        f3.close()

        # ANOVA
        otu_table_paths = [fp1,fp2]
        threshold = None
        _filter = 0
        otu_include = None
        otu_table1 = open(fp1,'U')
        otu_table2 = open(fp2,'U')
        category = 'cat1'

        # get expected ANOVA output from each file separately
        results1 = test_wrapper('ANOVA', otu_table1, category_mapping, category, threshold, \
                 _filter, otu_include=None)
        results2 = test_wrapper('ANOVA', otu_table2, category_mapping, category, threshold, \
                 _filter, otu_include=None)

        results = test_wrapper_multiple('ANOVA', otu_table_paths,
                                        category_mapping, category,
                                        threshold, _filter,
                                        otu_include)

        # multiple results should be combination of individual results
        results1 = [result.split('\t') for result in results1]        
        results1_dict = {}
        for line in results1[1:]:
            results1_dict[line[0]] = line[1:]
        results2 = [result.split('\t') for result in results2]
        results2_dict = {}
        for line in results2[1:]:
            results2_dict[line[0]] = line[1:]
        results = [result.split('\t') for result in results]
        results_dict = {}
        for line in results[1:]:
            results_dict[line[0]] = line[1:]
    
        # for every OTU in the combined results, compare to average of separate
        for k in results_dict.keys():
            # skip FDR corrected values because they are calculated _after_ the
            # method combines separate results
            for i in [0,1,3,4]:
                entry1 = float(results1_dict[k][i])
                entry2 = float(results2_dict[k][i])
                entry_combined = float(results_dict[k][i])
                mean = round((entry1+entry2)/2.0, 3)
                self.assertEqual(round(entry_combined, 3), mean)


        # correlation
        otu_table_paths = [fp1,fp2]
        threshold = None
        _filter = 0
        otu_include = None
        otu_table1 = open(fp1,'U')
        otu_table2 = open(fp2,'U')
        category = 'cat2'

        # get expected correlation output from each file separately
        results1 = test_wrapper('correlation', otu_table1, category_mapping, category, threshold, \
                 _filter, otu_include=None)
        results2 = test_wrapper('correlation', otu_table2, category_mapping, category, threshold, \
                 _filter, otu_include=None)

        results = test_wrapper_multiple('correlation', otu_table_paths,
                                        category_mapping, category,
                                        threshold, _filter,
                                        otu_include)

        # multiple results should be combination of individual results
        results1 = [result.split('\t') for result in results1]        
        results1_dict = {}
        for line in results1[1:]:
            results1_dict[line[0]] = line[1:]
        results2 = [result.split('\t') for result in results2]
        results2_dict = {}
        for line in results2[1:]:
            results2_dict[line[0]] = line[1:]
        results = [result.split('\t') for result in results]
        results_dict = {}
        for line in results[1:]:
            results_dict[line[0]] = line[1:]
    
        # for every OTU in the combined results, compare to average of separate
        for k in results_dict.keys():
            # skip FDR corrected values because they are calculated _after_ the
            # method combines separate results
            for i in [0,1,3,]:
                entry1 = float(results1_dict[k][i])
                entry2 = float(results2_dict[k][i])
                entry_combined = float(results_dict[k][i])
                mean = round((entry1+entry2)/2.0, 3)
                self.assertEqual(round(entry_combined, 3), mean)


        # G test
        otu_table_paths = [fp1,fp2]
        threshold = None
        _filter = 0
        otu_include = None
        otu_table1 = open(fp1,'U')
        otu_table2 = open(fp2,'U')
        category = 'cat1'

        # get expected G_TEST output from each file separately
        results1 = test_wrapper('g_test', otu_table1, category_mapping, category, threshold, \
                 _filter, otu_include=['0'])

        results2 = test_wrapper('g_test', otu_table2, category_mapping, category, threshold, \
                 _filter, otu_include=['0'])
        results = test_wrapper_multiple('g_test', otu_table_paths,
                                        category_mapping, category,
                                        threshold, _filter,
                                        otu_include=['0'])

        # multiple results should be combination of individual results
        results1 = [result.split('\t') for result in results1]        
        results1_dict = {}
        for line in results1[1:]:
            results1_dict[line[0]] = line[1:]
        results2 = [result.split('\t') for result in results2]
        results2_dict = {}
        for line in results2[1:]:
            results2_dict[line[0]] = line[1:]
        results = [result.split('\t') for result in results]
        results_dict = {}
        for line in results[1:]:
            results_dict[line[0]] = line[1:]
    
        # convert all numeric entries to floats rounded to 5 decimal places
        for k in results_dict.keys():
            # skip FDR corrected values because they are calculated _after_ the
            # method combines separate results
            for i in [0,1,2]:
                entry1 = float(results1_dict[k][i])
                entry2 = float(results2_dict[k][i])
                entry_combined = float(results_dict[k][i])
                mean = round((entry1+entry2)/2.0, 3)
                self.assertEqual(round(entry_combined, 3), mean)

        # clean up temporary files
        remove(fp1)
        remove(fp2)
        remove(fp3)


#run unit tests if run from command-line
if __name__ == '__main__':
    main()
