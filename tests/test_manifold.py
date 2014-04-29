#!/usr/bin/env python

__author__ = "Joshua Haas"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Joshua Haas,","Gregory Ditzler"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Joshua Haas"
__email__ = "laptopdude2@gmail.com"

from cogent.util.unit_test import TestCase, main
from qiime.manifold import compute_manifold

class ManifoldTests(TestCase):

    '''Tests of Top-Level Functions'''

    def setUp(self):
        self.biomstr =('{"id":null,"format": "Biological Observation Matrix 0.9.1-dev",' +
        '"format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",' +
        '"type": "OTU table","generated_by": "QIIME revision 1.4.0-dev",' +
        '"date": "2011-12-19T19:00:00","rows":[' +
        '{"id":"GG_OTU_1", "metadata":null},' +
        '{"id":"GG_OTU_2", "metadata":null},' +
        '{"id":"GG_OTU_3", "metadata":null},' +
        '{"id":"GG_OTU_4", "metadata":null},' +
        '{"id":"GG_OTU_5", "metadata":null},' +
        '{"id":"GG_OTU_6", "metadata":null}],"columns": [' +
        '{"id":"Sample1", "metadata":null},' +
        '{"id":"Sample2", "metadata":null},' +
        '{"id":"Sample3", "metadata":null},' +
        '{"id":"Sample4", "metadata":null},' +
        '{"id":"Sample5", "metadata":null},' +
        '{"id":"Sample6", "metadata":null}],' +
        '"matrix_type": "dense", "matrix_element_type": "int",' +
        '"shape": [6, 6],"data":[' +
        '[0,0,1,3,1,8],' +
        '[4,1,6,1,1,2],' +
        '[8,3,2,2,2,3],' +
        '[6,5,1,4,2,2],' +
        '[7,6,4,7,3,3],' +
        '[5,1,4,2,1,5]]}')

    def test_manifold(self):
        '''Manifold calculation should throw no errors'''
        result = compute_manifold([self.biomstr],"isomap")
        assert result

if __name__ == '__main__':
    main()
