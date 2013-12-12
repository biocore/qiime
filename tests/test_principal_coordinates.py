#!/usr/bin/env python


__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from qiime.principal_coordinates import pcoa
from cogent.util.unit_test import TestCase, main


class FunctionTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        self.distmtx_txt = """\tsam1\tsam2\tsam3
sam1\t0.00\t.18\t.44
sam2\t0.18\t0\t.66
sam3\t.44\t.66\t0""".split('\n')

    def test_pcoa(self):
        """ pcoa should throw no errors"""
        res = pcoa(self.distmtx_txt)
        assert res # formatting tested elsewhere


#run tests if called from command line
if __name__ == '__main__':
    main()