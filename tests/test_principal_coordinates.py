#!/usr/bin/env python


__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from unittest import TestCase, main
from StringIO import StringIO

from qiime.principal_coordinates import pcoa


class FunctionTests(TestCase):

    """Tests of top-level functions"""

    def setUp(self):
        self.distmtx_txt = StringIO(
        """\tsam1\tsam2\tsam3\n"""
        """sam1\t0.00\t0.18\t0.44\n"""
        """sam2\t0.18\t0.0\t.66\n"""
        """sam3\t0.44\t0.66\t0.0""")

    def test_pcoa(self):
        """ pcoa should throw no errors"""
        res = pcoa(self.distmtx_txt)
        assert res  # formatting tested elsewhere


# run tests if called from command line
if __name__ == '__main__':
    main()
