#!/usr/bin/env python
# File created on 19 May 2011


from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren, Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"
__status__ = "Release"
 


from cogent.util.unit_test import TestCase,main
from qiime.parse import parse_mapping_file_to_dict, parse_rarefaction
from numpy import array
from qiime.compare_alpha_diversity import (compare_alpha_diversities,
    extract_rarefaction_scores_at_depth,
    make_SampleIds_rarefaction_columns_dict,
    make_category_values_Id_dict,
    make_value_pairs_from_category,
    map_category_value_pairs_to_Ids,
    convert_SampleIds_to_rarefaction_mtx)

class TopLevelTests(TestCase):
    """Tests of top level functions"""
    
    def setUp(self):
        """define data for tests"""
        self.rarefaction_file = \
         ['\tsequences per sample\titeration\t123\t234\t345',
          'rare10.txt\t10\t0\t1.99181\t5.42877\t2.13996',
          'rare10.txt\t10\t1\t2.07163\t1.42877\t2.37055',
          'rare310.txt\t310\t0\t8.83115\t6.42877\t11.00725',
          'rare310.txt\t310\t1\t10.05242\t9.42877\t8.24474',
          'rare810.txt\t810\t0\t12.03067\tn/a\t11.58928',
          'rare910.txt\t910\t1\t12.9862\t2.42877\t11.58642']
        
        self.rarefaction_data = parse_rarefaction(self.rarefaction_file)
        
        self.mapping_file = \
        ['#SampleID\tTreatment\tLinker'+\
         'PrimerSequence\tDose\tTTD\tDescription',
         '#Comment Line',
         '123\tAAAA\tBBBB\tHigh\t31\tM_ID_123',
         '234\tCCCC\tDDDD\tLow\t67\tM_ID_234',
         '345\tAAAA\tFFFF\tMed\t21\tM_ID_345']
        
        self.mapping_data = \
         parse_mapping_file_to_dict(self.mapping_file)[0]
        self.value_pairs_Dose = \
         [('Low','Med'),('Low','High'),('Med','High')]                         
        self.value_pairs_TTD = \
         [('67', '21'), ('67', '31'), ('21', '31')]
        self.value_pairs_Treatment = \
         [('CCCC', 'AAAA')]
        self.cat_val_Dose = \
         {'High': ['123'], 'Low': ['234'], 'Med': ['345']}
        self.cat_val_TTD = \
         {'21': ['345'], '31': ['123'], '67': ['234']}
        self.cat_val_Treatment = \
         {'AAAA': ['345', '123'], 'CCCC': ['234']}
        self.Id_pairs_Dose = \
         [(['234'], ['345']), (['234'], ['123']), (['345'], ['123'])]
        self.Id_pairs_TTD = \
         [(['234'], ['345']), (['234'], ['123']), (['345'], ['123'])]
        
        self.Id_pairs_Treatment = \
         [(['234'], ['345', '123'])]
        
        self.rarefaction_cols_dict = \
         {'123': 0, '234': 1, '345': 2}
       
        self.extracted_mtx_10 = \
         array([[ 1.99181,  5.42877,  2.13996],
               [ 2.07163,  1.42877,  2.37055]])
        
        self.extracted_mtx_310 = \
         array([[  8.83115,   6.42877,  11.00725],
               [ 10.05242,   9.42877,   8.24474]])
        
        self.extracted_mtx_910 = \
         array([[ 12.9862 ,   2.42877,  11.58642]])
        
        self.sample_pair1 = \
         (['234'], ['345', '123'])
            
        self.rarefaction_mtx_for_sample_pair1_0 = \
         array([[ 5.42877],
               [ 1.42877]])
        
        self.rarefaction_mtx_for_sample_pair1_1 = \
         array([[ 2.13996,  1.99181],
               [ 2.37055,  2.07163]])
        
        self.compared_alpha_diversities_TTD = \
         {'TTD':{('21','31'):(1.8321466933860993,0.20839398129924847),
         ('67', '21'): (0.58578495700890432, 0.61731739324369639),
         ('67', '31'): (0.69838596448703294, 0.55721515283248324)}}
        
    def test_make_value_pairs_from_category(self):
        """check value pairs returns correct unique pairs for categories
        """
        
          
        self.assertEqual(
         self.value_pairs_Dose,
         make_value_pairs_from_category(self.mapping_data, 'Dose'))
        
        self.assertEqual(
            self.value_pairs_TTD, 
            make_value_pairs_from_category(self.mapping_data, 'TTD'))
        
        self.assertEqual(
            self.value_pairs_Treatment, 
            make_value_pairs_from_category(self.mapping_data, 
                                            'Treatment'))
        
        self.assertRaises(
         ValueError,
         make_value_pairs_from_category, self.mapping_data, 'WrongCat')
                         
    def test_make_category_values_Id_dict(self):
        """check value pairs reference correct Id pairs"""
             
        self.assertEqual(
            self.cat_val_Dose,
            make_category_values_Id_dict(self.mapping_data, 'Dose'))
        
        self.assertEqual(
            self.cat_val_TTD,
            make_category_values_Id_dict(self.mapping_data, 'TTD'))
        
        self.assertEqual(
            self.cat_val_Treatment,
            make_category_values_Id_dict(self.mapping_data,
                                         'Treatment'))
    
    def test_map_category_value_pairs_to_Ids(self):
        """check value pairs converted to correct Ids"""
                
        self.assertEqual(
            self.Id_pairs_Dose,
            map_category_value_pairs_to_Ids(self.value_pairs_Dose,
                                            self.cat_val_Dose))

        self.assertEqual(
            self.Id_pairs_TTD,
            map_category_value_pairs_to_Ids(self.value_pairs_TTD,
                                            self.cat_val_TTD))
        
        self.assertEqual(
            self.Id_pairs_Treatment,
            map_category_value_pairs_to_Ids(
            self.value_pairs_Treatment,
            self.cat_val_Treatment))
    
    def test_make_SampleIds_rarefaction_columns_dict(self):
        """ """
                
        self.assertEqual(
         self.rarefaction_cols_dict,
         make_SampleIds_rarefaction_columns_dict(self.rarefaction_data))
        
    def test_extract_rarefaction_scores_at_depth(self):
        """check correct errors raised for wrong depths and correct mtx
        """
           
        self.assertEqual(
         self.extracted_mtx_10,
         extract_rarefaction_scores_at_depth(10,self.rarefaction_data))
        
        self.assertEqual(
         self.extracted_mtx_310,
         extract_rarefaction_scores_at_depth(310,self.rarefaction_data))
        
        self.assertEqual(
         self.extracted_mtx_910,
         extract_rarefaction_scores_at_depth(910,self.rarefaction_data))
        
        self.assertRaises(
         ValueError,
         extract_rarefaction_scores_at_depth,810,self.rarefaction_data)
        
        self.assertRaises(
         ValueError,
         extract_rarefaction_scores_at_depth,100,self.rarefaction_data)
        
    def test_convert_SampleIds_to_rarefaction_mtx(self):
        """check correct reduced rarefaction scores mtx produced"""
                 
        self.rarefaction_mtx_for_sample_pair1_0 = \
         array([[ 5.42877],
               [ 1.42877]])
        
        self.rarefaction_mtx_for_sample_pair1_1 = \
         array([[ 2.13996,  1.99181],
               [ 2.37055,  2.07163]])
        
        self.assertEqual(
            convert_SampleIds_to_rarefaction_mtx(
                self.sample_pair1[0],
                self.extracted_mtx_10,
                self.rarefaction_cols_dict),
            self.rarefaction_mtx_for_sample_pair1_0)
        
        self.assertEqual(
            convert_SampleIds_to_rarefaction_mtx(
                self.sample_pair1[1],
                self.extracted_mtx_10,
                self.rarefaction_cols_dict),
            self.rarefaction_mtx_for_sample_pair1_1)
        
    def test_compare_alpha_diversity(self):
        """test main function properly compares alpha diversities"""
                
        self.assertEqual(
            self.compared_alpha_diversities_TTD,
            compare_alpha_diversities(self.rarefaction_file,
                                      self.mapping_file, 'TTD', 10))


if __name__ == "__main__":
    main()
