#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division

__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler", "Calvin Morrison"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Release"
 

import qiime.fizzy as fizzy 

def main():
	"""
		main()
		No input agruments

		Run the fizzy module with data that is in qiime/tests/test_fizzy/
		directory. The result should be in output.txt. If the file is not
		there then there was likely an error. 
	"""
	data_file = open('./test_fizzy/data.biom','U')
	label_file = open('./test_fizzy/labels.tsv','U')
	n_select = 25
	method = 'mim'
	outfile = './test_fizzy/output.txt'
	try:
		fizzy.run_feature_selection(data_file, label_file, outfile, method, n_select)
		print 'Refer to the test_fizzy directory for output.txt '
	except ValueError:
		print 'Something is broken in Fizzy!'

	return None


 if __name__ == "__main__":
 	main()
