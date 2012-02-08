#!/usr/bin/env python
# File created on 08 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from qiime.quality import ascii_to_phred33, ascii_to_phred64, ascii_to_phred

class QualityTests(TestCase):
    
    def setUp(self):
        """ """
        pass
    
    def test_ascii_to_phred33(self):
        """ conversions from ascii to phred function as expected """
        self.assertEqual(ascii_to_phred33('!'),0)
        self.assertEqual(ascii_to_phred33('?'),30)

    def test_ascii_to_phred64(self):
        """ conversions from ascii to phred function as expected """
        self.assertEqual(ascii_to_phred64('@'),0)
        self.assertEqual(ascii_to_phred64('^'),30)

    def test_ascii_to_phred(self):
        """ conversions from ascii to phred function as expected """
        self.assertEqual(ascii_to_phred('x',120),0)
        self.assertEqual(ascii_to_phred('x',119),1)


if __name__ == "__main__":
    main()