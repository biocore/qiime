#!/usr/bin/env python
# File created on 27 Sep 2011
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena", "Andrew J. King", "Michael S. Robeson",]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"

from qiime.distance_matrix_from_mapping import compute_distance_matrix_from_metadata, dist_vincenty, calculate_dist_vincenty
from numpy import array
from cogent.util.unit_test import TestCase, main
import StringIO
    

class FunctionTests(TestCase):
  """Tests of top-level functions"""

  def setUp(self):
    self.fasting_map = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Float_Col	Description
#Example mapping file for the QIIME analysis package. These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	.1	Control_mouse__I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	.2	Control_mouse__I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	20061126	.3	Control_mouse__I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	.4	Control_mouse__I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	.5	Control_mouse__I.D._593
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	.6	Fasting_mouse__I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	.7	Fasting_mouse__I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	.8	Fasting_mouse__I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	.9	Fasting_mouse__I.D._636"""
    self.DOB = [20061218, 20061218, 20061126, 20070314, 20071210, 20071112, 20080116, 20080116, 20080116]
    self.Float_Col = [.1, .2, .3, .4, .5, .6, .7, .8, .9]
    self.latitudes = [30, 20, 30, 30, 0, 1, 90, 89, 0, 0]
    self.longitudes = [60, -50, 60, 60, 0, 0, 0, 0, 0, 0]

  def test_compute_distance_matrix_from_metadata_int(self):
    """ distance calculations on ints should throw no errors"""
    exp_out = array([[0, 0, 92, 9096, 9992, 9894, 18898, 18898, 18898], [0, 0, 92, 9096, 9992, 9894, 18898, 18898, 18898],
      [92, 92, 0, 9188, 10084, 9986, 18990, 18990, 18990], [9096, 9096, 9188, 0, 896, 798, 9802, 9802, 9802],
      [9992, 9992, 10084, 896, 0, 98, 8906, 8906, 8906], [9894, 9894, 9986, 798, 98, 0, 9004, 9004, 9004],
      [18898, 18898, 18990, 9802, 8906, 9004, 0, 0, 0], [18898, 18898, 18990, 9802, 8906, 9004, 0, 0, 0],
      [18898, 18898, 18990, 9802, 8906, 9004, 0, 0, 0]])

    res_out = compute_distance_matrix_from_metadata(self.DOB)
    self.assertFloatEqual(exp_out, res_out)
    
  def test_compute_distance_matrix_from_metadata_floats(self):
    """ distance calculations on floats should throw no errors"""
    # testing floats
    exp_out = array([[0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], [0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
      [0.2, 0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6], [0.3, 0.2, 0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5],
      [0.4, 0.3, 0.2, 0.1, 0., 0.1, 0.2, 0.3, 0.4], [0.5, 0.4, 0.3, 0.2, 0.1, 0., 0.1, 0.2, 0.3],
      [0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0., 0.1, 0.2], [0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0., 0.1],
      [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.]])
 
    res_out = compute_distance_matrix_from_metadata(self.Float_Col)
    self.assertFloatEqual(exp_out, res_out)
   
  def test_dist_vincenty(self):
    """dist_Vincenty:Returns distance in meters between two lat long points"""
    lat1, lon1, lat2, lon2, expected_value = 30, 60, 20, -50, 10709578.387
    value = dist_vincenty(lat1, lon1, lat2, lon2, 20)
    self.assertFloatEqual(value, expected_value)
    
    lat1, lon1, lat2, lon2, expected_value = 30, 60, 30, 60, 0
    value = dist_vincenty(lat1, lon1, lat2, lon2, 20)
    self.assertFloatEqual(value, expected_value)
    
    lat1, lon1, lat2, lon2, expected_value = 0,  0,  1, 0, 110574.389
    value = dist_vincenty(lat1, lon1, lat2, lon2, 20)
    self.assertFloatEqual(value, expected_value)
    
    lat1, lon1, lat2, lon2, expected_value = 90,  0, 89, 0, 111693.865
    value = dist_vincenty(lat1, lon1, lat2, lon2, 20)
    self.assertFloatEqual(value, expected_value)
    
    lat1, lon1, lat2, lon2, expected_value = 90,  0, -90, 0, 20003931.459
    value = dist_vincenty(lat1, lon1, lat2, lon2, 20)
    self.assertFloatEqual(value, expected_value)
    
    lat1, lon1, lat2, lon2, expected_value = 90,  0,  0, 0, 10001965.729
    value = dist_vincenty(lat1, lon1, lat2, lon2, 20)
    self.assertFloatEqual(value, expected_value)
    
    lat1, lon1, lat2, lon2, expected_value = 0,  0,  0, 0, 0
    value = dist_vincenty(lat1, lon1, lat2, lon2, 20)
    self.assertFloatEqual(value, expected_value)
    
    # test for not converge
    lat1, lon1, lat2, lon2 = 0, 180,  0, 0
    self.assertRaises(ValueError, dist_vincenty, lat1, lon1, lat2, lon2, 20)

  def test_calculate_dist_vincenty(self):
    exp_out = array([[0.0, 10709578.387, 0.0, 0.0, 7154900.607, 7094106.828, 6681852.331, 6626434.332, 7154900.607, 7154900.607],
        [10709578.387, 0.0, 10709578.387, 10709578.387, 5877643.846, 5831009.412, 7789599.475, 7718017.604, 5877643.846, 5877643.846],
        [0.0, 10709578.387, 0.0, 0.0, 7154900.607, 7094106.828, 6681852.331, 6626434.332, 7154900.607, 7154900.607],
        [0.0, 10709578.387, 0.0, 0.0, 7154900.607, 7094106.828, 6681852.331, 6626434.332, 7154900.607, 7154900.607],
        [7154900.607, 5877643.846, 7154900.607, 7154900.607, 0.0, 110574.389, 10001965.729, 9890271.864, 0.0, 0.0],
        [7094106.828, 5831009.412, 7094106.828, 7094106.828, 110574.389, 0.0, 9891391.341, 9779697.476, 110574.389, 110574.389],
        [6681852.331, 7789599.475, 6681852.331, 6681852.331, 10001965.729, 9891391.341, 0.0, 111693.865, 10001965.729, 10001965.729],
        [6626434.332, 7718017.604, 6626434.332, 6626434.332, 9890271.864, 9779697.476, 111693.865, 0.0, 9890271.864, 9890271.864],
        [7154900.607, 5877643.846, 7154900.607, 7154900.607, 0.0, 110574.389, 10001965.729, 9890271.864, 0.0, 0.0],
        [7154900.607, 5877643.846, 7154900.607, 7154900.607, 0.0, 110574.389, 10001965.729, 9890271.864, 0.0, 0.0]])
    
    res_out = calculate_dist_vincenty(self.latitudes, self.longitudes)
    
    self.assertFloatEqual(res_out, exp_out)
    

#run tests if called from command line
if __name__ == '__main__':
  main()
  
  
  

