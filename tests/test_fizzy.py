#!/usr/bin/env python
from __future__ import division
import qiime.fizzy as fizzy 
from cogent.util.unit_test import TestCase, main
import StringIO
import numpy as np

__author__ = ["Gregory Ditzler", "Calvin Morrision"]
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler", "Calvin Morrison", "Gail Rosen"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Development"

class FizzyTests(TestCase):

    def setUp(self):
        self.biom_file_handle = StringIO.StringIO(biom_file_string)
        self.map_file_handle = StringIO.StringIO(map_file_string)

    def uniform_data(self, n_observations, n_features, n_select):
        """ Generate some uniform data to use for test_pyfeast_run """
        xmax = 10
        xmin = 1

        data = 1.0*np.random.randint(xmax + 1, size = (n_features, n_observations))
        labels = np.zeros(n_observations)
        delta = n_select * (xmax - xmin) / 2.0

        for m in range(n_observations):
            zz = 0.0
            for k in range(n_select):
                zz += data[k, m]
                if zz > delta:
                    labels[m] = 1
                else:
                    labels[m] = 2
        data = data.transpose()
        return data, labels


    def test_get_fs_methods(self):
        """this test is going to make sure that we are only implementing feature selection methods that are in a predefined list. 
        """
        self.assertEqual(fizzy.get_fs_methods(), information_theoretic_methods)

    def test_parse_biom(self):
        """this test is going to ensure that we can properly parse our biom file"""
        correct_biom = ([[  1.,   0.,   1.,   6.], [  0.,   5.,   1.,  10.], [  4.,   7.,   9.,   8.]], [u'OTU0', u'OTU1', u'OTU2', u'OTU3'], [u'ID0', u'ID1', u'ID2'])

        parsed_biom = fizzy.parse_biom(self.biom_file_handle)
        self.assertEqual(parsed_biom, correct_biom)
            
    
    def test_parse_map_file(self):
        """this test ensures that the map file is properly parsed"""

        correct_map = [0,1.,1.]
        parsed_map = fizzy.parse_map_file(self.map_file_handle, "Class", [u'ID0',u'ID1',u'ID2'])
        self.assertEqual(parsed_map, correct_map)


    def test_parse_map_invalid_column_name(self):
        """test the functionality of the map parser."""
        self.assertRaises(ValueError, fizzy.parse_map_file, self.map_file_handle, "Unknown", [u'ID0',u'ID1',u'ID2'])

    def test_parse_map_invalid_sample_name(self):
        """test with invalid sample names"""
        self.assertRaises(ValueError, fizzy.parse_map_file, self.map_file_handle, "Class", [u'ID99',u'ID1',u'ID2'])
    
    def test_parse_map_invalid_equal_num_sample_classes(self):
        """test with an invalid number of samples"""
        self.assertRaises(ValueError, fizzy.parse_map_file, self.map_file_handle, "Class", [u'ID0',u'ID1'])

    def test_run_pyfeast(self):
        """assert that each of the selected features are in the appropriate range"""
        n_select = 5
        (data, labels)  = self.uniform_data(500, 50, n_select)
        variable_names = range(len(labels)) # we don't need this to  be differenttest pyfeast, so just give it the same sized vector

        selected_features = fizzy.run_pyfeast(data, labels, variable_names, "MIM", n_select)
        for k in selected_features:
          self.assertTrue(k in range(n_select))

    def test_run_pyfeast_invalid(self):
        """test the run_pyfeast function with an invalid algorithm name"""

        self.assertRaises(AttributeError, fizzy.run_pyfeast, None, None, None, "An Invalid Name", None)

    def test_run_feature_selection(self):
        """test the feature selection on a toy problem"""
        column_name = "Class"

        selected_features = fizzy.run_feature_selection(self.biom_file_handle, self.map_file_handle, column_name, n_select=3)
        self.assertEqual(selected_features, ['OTU0', "OTU1", "OTU3"])

    def test_zero_selection(self):
        """assert that fizzy throws an error for zero features being selected"""
        (data, labels)  = self.uniform_data(500, 50, n_select)
        variable_names = range(len(labels)) # we don't need this to  be differenttest pyfeast, so just give it the same sized vector
        self.assertRaises(ValueError, fizzy.run_pyfeast, data, labels, variable_names, "MIM", 0)

    def test_negative_selection(self):
        """assert that fizzy throws an error for negative features being selected"""
        (data, labels)  = self.uniform_data(500, 50, n_select)
        variable_names = range(len(labels)) # we don't need this to  be differenttest pyfeast, so just give it the same sized vector
        self.assertRaises(ValueError, fizzy.run_pyfeast, data, labels, variable_names, "MIM", -1)
     
    def test_too_large_selection(self):
        """assert that fizzy throws an error for requesting more features than exist."""
        (data, labels)  = self.uniform_data(500, 50, n_select)
        variable_names = range(len(labels)) # we don't need this to  be differenttest pyfeast, so just give it the same sized vector
        slef.assertRaises(ValueError, fizzy.run_pyfeast, data, labels, variable_names, "MIM", 10000)
     
information_theoretic_methods = ['CIFE','CMIM','CondMI', 'Condred','ICAP','JMI','MIM','MIFS','mRMR']
map_file_string = """#SampleID\tClass\nID0\tS1\nID1\tS2\nID2\tS2\n"""
biom_file_string = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date":"2013-03-27T13:59:38.949014","matrix_type": "sparse","matrix_element_type":"float","shape": [4, 3],     "data": [[0,0,1.0],[0,2,4.0],[1,1,5.0],[1,2,7.0],[2,0,1.0],[2,1,1.0],[2,2,9.0],[3,0,6.0],[3,1,10.0],[3,2,8.0]],"rows":[{"id": "OTU0", "metadata": null},{"id": "OTU1", "metadata": null},{"id": "OTU2", "metadata": null},{"id": "OTU3", "metadata": null}],"columns": [{"id": "ID0", "metadata": null},{"id": "ID1", "metadata": null},{"id": "ID2", "metadata": null}]}"""

if __name__ == "__main__":
    main()
