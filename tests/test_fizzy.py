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
		"""
		# it should be safe to assume the user is working with 
		# linux / unix. this will not be supported by windows. 

		fizzy.test_parse_biom()
		return None 
	
	def test_parse_map_file(self):
		"""
		"""
		return None 

	def test_run_pyfeast():
		"""
		"""
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
		return None 

		
information_theoretic_methods = ['CIFE','CMIM','CondMI', 'Condred','ICAP','JMI','MIM','MIFS','mRMR']
biome_file = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-03-27T13:59:38.949014","matrix_type": "sparse","matrix_element_type": "float","shape": [4, 3],     "data": [[0,0,1.0],[0,2,4.0],[1,1,5.0],[1,2,7.0],[2,0,1.0],[2,1,1.0],[2,2,9.0],[3,0,6.0],[3,1,10.0],[3,2,8.0]],"rows": [{"id": "OTU0", "metadata": null},{"id": "OTU1", "metadata": null},{"id": "OTU2", "metadata": null},{"id": "OTU3", "metadata": null}],"columns": [{"id": "S1", "metadata": null},{"id": "S2", "metadata": null},{"id": "S3", "metadata": null}]}"""


if __name__ == "__main__":
	main()
