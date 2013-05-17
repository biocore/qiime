#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division


"""
	This is a test unit for the Fizzy feature selection algorithm, which uses
	the PyFeast python module. If the unit test fails, then it is likely that 
	you do not have PyFeast installed on the system. 

	After PyFeast is installed you should be good to go. If you have any 
	questions feel free to email gregory.ditzler@gmail.com 
"""

"""
test_fizzy.py
Copyright (C) 2013 Gregory Ditzler

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler", "Calvin Morrison", "Gail Rosen"]
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
	label_file = open('./test_fizzy/labels.map','U')
	n_select = 25
	method = 'mim'
	col_name = 'Class'
	outfile = './test_fizzy/output.txt'
	try:
		# if this throws an error then it is likely that PyFeast has not
		# been installed. 
		fizzy.run_feature_selection(data_file, label_file, col_name, outfile, method, n_select)
		print 'Refer to the test_fizzy directory for output.txt '
		f = open('./test_fizzy/output.txt','U')
		print f.read()
		f.close()
	except ValueError:
		print 'Something is broken in Fizzy! Do you have PyFeast installed?'

	return None


if __name__ == "__main__":
	main()
