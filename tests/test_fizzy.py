#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division

__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler", "Calvin Morrison", "Gail Rosen"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Development"

import qiime.fizzy as fizzy 
import sys, os 
from cogent.util.unit_test import TestCase, main
from qiime.util import (get_qiime_scripts_dir,
	load_qiime_config,
	get_tmp_filename,
	get_qiime_temp_dir,
	create_dir)
import shutil
import StringIO

class TestFizzy(TestCase):
	"""
		docstring for TestFizzy
	"""

	def setUp(self):
		"""
			setup function defined by the bioler plate qiime 
			test functions. 
		"""
		self.dirs_to_remove = []
		self.files_to_remove = []

		tmp_dir = get_qiime_temp_dir()
		self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
			prefix='qiime_parallel_tests_',
			suffix='',
			result_constructor=str)
		self.dirs_to_remove.append(self.test_out)
		create_dir(self.test_out)

	def tearDown(self):
		"""
			tear down function defined by the bioler plate qiime 
			test functions. 
		"""
		for d in self.dirs_to_remove:
			if os.path.exists(d):
				shutil.rmtree(d)
	
	def test_get_fs_methods(self):
		"""
			this test is going to make sure that we are only 
			implementing feature selection methods that are in 
			a predefined list. 
		"""
		self.assertEqual(information_theoretic_methods, 
			fizzy.get_fs_methods())
		return None 

	def test_parse_biom(self):
		"""
      this test is going to ensure that we can properly parse our biom file
		"""
		# it should be safe to assume the user is working with 
		# linux / unix. this will not be supported by windows. 

		correct_biom = ([[  1.,   0.,   1.,   6.], [  0.,   5.,   1.,  10.], [  4.,   7.,   9.,   8.]], [u'OTU0', u'OTU1', u'OTU2', u'OTU3'], [u'ID0', u'ID1', u'ID2'])

		biom_file_handle = StringIO.StringIO(biom_file_string)
		parsed_biom = fizzy.parse_biom(biom_file_handle)
		self.assertEqual(correct_biom, parsed_biom)
			
		return None 
	
	def test_parse_map_file(self):
		"""
			this test ensures that the map file is proerply parsed
		"""

		correct_map = [0,1.,1.]

		map_file_handle = StringIO.StringIO(map_file_string)
		parsed_map = fizzy.parse_map_file(map_file_handle, "Class", [u'ID0',u'ID1',u'ID2'])

		self.assertEqual(correct_map, parsed_map)
		return None 

	def uniform_data(self, n_observations, n_features, n_select):
		""" Generate some uniform data to use for test_pyfeast_run """

		import numpy as np
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

	def test_run_pyfeast(self):
		"""
			test the run_pyfeast function by creating some uniform data, and running
			the MIM feature selection method.

			Assert that each of the selected features are in the appropriate range
		"""
		n_select = 5
		(data, labels)  = self.uniform_data(100, 50, n_select)
		variable_names = labels # we don't need this to  be differenttest pyfeast, so just give it the same sized vector
		
		selected_features = fizzy.run_pyfeast(data, labels, variable_names, "MIM", n_select)
		selected_features = sorted(selected_features)

		for k in range(n_select):
			if k != selected_features[k]:
				self.assertEqual(True, k in range(n_select))

		return None 

	def test_write_output_file(self):
		"""
			wrtie a temporary file to 
		"""
		fizzy.write_output_file("test",self.test_out + "/test.txt")
		self.assertEqual(open(self.test_out + "/test.txt").read().replace("\n","").split(),
			["test"])
		return None 

	def test_run_feature_selection(self):
		"""
		"""
  
		map_file_handle = StringIO.StringIO(map_file_string)
		biom_file_handle = StringIO.StringIO(biom_file_string)
		column_name = "Class"
		method = "MIM"

		fizzy.run_feature_selection(biom_file_handle, map_file_handle, column_name, self.test_out + "fs_test", method='MIM', n_select=15)

		return None 

		
information_theoretic_methods = ['CIFE','CMIM','CondMI', 'Condred','ICAP','JMI','MIM','MIFS','mRMR']
map_file_string = """#SampleID	Class\nID0	S1\nID1	S2\nID2	S2\n"""
biom_file_string = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date":"2013-03-27T13:59:38.949014","matrix_type": "sparse","matrix_element_type":"float","shape": [4, 3],     "data": [[0,0,1.0],[0,2,4.0],[1,1,5.0],[1,2,7.0],[2,0,1.0],[2,1,1.0],[2,2,9.0],[3,0,6.0],[3,1,10.0],[3,2,8.0]],"rows":[{"id": "OTU0", "metadata": null},{"id": "OTU1", "metadata": null},{"id": "OTU2", "metadata": null},{"id": "OTU3", "metadata": null}],"columns": [{"id": "ID0", "metadata": null},{"id": "ID1", "metadata": null},{"id": "ID2", "metadata": null}]}"""

if __name__ == "__main__":
	main()
