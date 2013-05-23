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


# def main():
# 	"""
# 		main()
# 		No input agruments

# 		Run the fizzy module with data that is in qiime/tests/test_fizzy/
# 		directory. The result should be in output.txt. If the file is not
# 		there then there was likely an error. 
# 	"""
# 	data_file = open('../qiime_test_data/fizzy/data.biom','U')
# 	label_file = open('../qiime_test_data/fizzy/labels.map','U')
# 	n_select = 25
# 	method = 'MIM'
# 	col_name = 'Class'
# 	outfile = '../qiime_test_data/fizzy/output.txt'
# 	try:
# 		# if this throws an error then it is likely that PyFeast has not
# 		# been installed. 
# 		fizzy.run_feature_selection(data_file, label_file, col_name, outfile, method, n_select)
# 		print 'Refer to the test_fizzy directory for output.txt '
# 		f = open('./test_fizzy/output.txt','U')
# 		print f.read()
# 		f.close()
# 	except:
# 		raise RuntimeError("""Error :: tests/test_fizzy.py : There was an
# 		error running the fizzy test script. It is likely that the PyFeast
# 		Python module is not installed on the system. """)
		

# 	return None

class TestFizzy(TestCase):
	"""docstring for TestFizzy"""

	def setUp(self):
		"""
			set up
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
			tear donw
		"""
		for d in self.dirs_to_remove:
			if os.path.exists(d):
				shutil.rmtree(d)
	
	def test_script(self):
		"""
			test fizzy 
		"""

		

if __name__ == "__main__":
	main()
