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
	"""docstring for TestFizzy"""

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
	
	def test_get_fs_methods():
		"""
			this test is going to make sure that we are only 
			implementing feature selection methods that are in 
			a predefined list. 
		"""
		information_theoretic_methods = ['CIFE','CMIM','CondMI',\
			'Condred','ICAP','JMI','MIM','MIFS','mRMR']
		self.assertEqual(information_theoretic_methods, 
			fizzy.get_fs_methods())
		return None 

	def test_parse_biom():
		"""
		"""
		return None 
	
	def test_parse_map_file():
		"""
		"""
		return None 

	def test_run_pyfeast():
		"""
		"""
		return None 

	def test_write_output_file():
		"""
		"""
		return None 

	def test_run_feature_selection():
		"""
		"""
		return None 

		

if __name__ == "__main__":
	main()
